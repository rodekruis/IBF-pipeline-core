from __future__ import annotations

import os.path
import copy
import time
import os
import json
import cdsapi
from ibfpipeline.core.secrets import Secrets
from ibfpipeline.core.settings import Settings
from ibfpipeline.drought.data import (
    AdminDataSet,
    AdminDataUnit,
    ForecastDataUnit,
    ClimateRegionDataSet,
    ClimateRegionDataUnit,
)
from urllib.error import HTTPError
import urllib.request, json
from datetime import datetime, timedelta, date
import azure.cosmos.cosmos_client as cosmos_client
import logging
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import requests
import geopandas as gpd
import shutil
from itertools import product
from azure.storage.blob import BlobServiceClient
from azure.core.exceptions import ResourceNotFoundError

COSMOS_DATA_TYPES = [
    "climate-region",
    "seasonal-rainfall-forecast",
    "seasonal-rainfall-forecast-climate-region",
]


def get_cosmos_query(
    start_date=None,
    end_date=None,
    country=None,
    adm_level=None,
    climate_region_code=None,
    pcode=None,
    lead_time=None,
):
    query = "SELECT * FROM c WHERE "
    if start_date is not None:
        query += f'c.timestamp >= "{start_date.strftime("%Y-%m-%dT%H:%M:%S")}" '
    if end_date is not None:
        query += f'AND c.timestamp <= "{end_date.strftime("%Y-%m-%dT%H:%M:%S")}" '
    if country is not None:
        query += f'AND c.country = "{country}" '
    if adm_level is not None:
        query += f'AND c.adm_level = "{adm_level}" '
    if climate_region_code is not None:
        query += f'AND c.climate_region_code = "{climate_region_code}" '
    if pcode is not None:
        query += f'AND c.adm_level = "{pcode}" '
    if lead_time is not None:
        query += f'AND c.adm_level = "{lead_time}" '
    if query.endswith("WHERE "):
        query = query.replace("WHERE ", "")
    query = query.replace("WHERE AND", "WHERE")
    return query


def get_data_unit_id(data_unit: AdminDataUnit, dataset: AdminDataSet):
    """Get data unit ID"""
    if hasattr(data_unit, "pcode") and getattr(data_unit, "pcode") is not None:
        if hasattr(data_unit, "lead_time"):
            id_ = f"{data_unit.pcode}_{dataset.timestamp.strftime('%Y-%m-%dT%H:%M:%S')}_{data_unit.lead_time}"
        else:
            id_ = f"{data_unit.pcode}_{dataset.timestamp.strftime('%Y-%m-%dT%H:%M:%S')}"
    elif hasattr(data_unit, "climate_region_code"):
        if hasattr(data_unit, "lead_time"):
            id_ = f"{data_unit.climate_region_code}_{dataset.timestamp.strftime('%Y-%m-%dT%H:%M:%S')}_{data_unit.lead_time}"
        else:
            id_ = f"{data_unit.climate_region_code}_{dataset.timestamp.strftime('%Y-%m-%dT%H:%M:%S')}"
    else:
        id_ = f"{dataset.timestamp.strftime('%Y-%m-%dT%H:%M:%S')}"
    return id_


def forecast_trigger_status(triggered: bool, trigger_class: str):
    """determine if forecast is a trigger for IBF portal if trigger status is true and trigger activation is enabled in config file the
    trigger staus will be 1 , else 0"""
    if triggered:
        if trigger_class == "enabled":
            return 1
        else:
            return 0
    else:
        return 0


class Load:
    """Download/upload data from/to a data storage"""

    def __init__(self, settings: Settings = None, secrets: Secrets = None):
        self.secrets = None
        self.settings = None
        if settings is not None:
            self.set_settings(settings)
        if secrets is not None:
            self.set_secrets(secrets)
        self.rasters_sent = []

    def set_settings(self, settings):
        """Set settings"""
        if not isinstance(settings, Settings):
            raise TypeError(f"invalid format of settings, use settings.Settings")
        settings.check_settings(
            ["postgresql_server", "postgresql_port", "postgresql_database"]
        )
        self.settings = settings

    def set_secrets(self, secrets):
        """Set secrets for storage"""
        if not isinstance(secrets, Secrets):
            raise TypeError(f"invalid format of secrets, use secrets.Secrets")
        secrets.check_secrets(
            [
                "COSMOS_URL",
                "COSMOS_KEY",
                "BLOB_ACCOUNT_NAME",
                "BLOB_ACCOUNT_KEY",
                "IBF_API_URL",
                "IBF_API_USER",
                "IBF_API_PASSWORD",
                "CDSAPI_KEY",
            ]
        )
        self.secrets = secrets

    def get_population_density(self, country: str, file_path: str):
        """Get population density data from worldpop and save to file_path"""
        r = requests.get(
            f"{self.settings.get_setting('worldpop_url')}/{country.upper()}/{country.lower()}_ppp_2020_UNadj_constrained.tif"  # f"{self.settings.get_setting('worldpop_url')}/{country.upper()}/{country.lower()}_ppp_2022_1km_UNadj_constrained.tif"
        )
        if "404 Not Found" in str(r.content):
            raise FileNotFoundError(
                f"Population density data not found for country {country}"
            )
        with open(file_path, "wb") as file:
            file.write(r.content)

    def get_adm_boundaries(self, country: str, adm_level: int) -> gpd.GeoDataFrame:
        """Get admin areas from IBF API"""
        try:
            adm_boundaries = self.ibf_api_get_request(
                f"admin-areas/{country}/{adm_level}",
            )
            gdf_adm_boundaries = gpd.GeoDataFrame.from_features(
                adm_boundaries["features"]
            )
            gdf_adm_boundaries.set_crs(epsg=4326, inplace=True)
        except HTTPError:
            raise FileNotFoundError(
                f"Admin areas for country {country}"
                f" and admin level {adm_level} not found"
            )
        return gdf_adm_boundaries

    def __ibf_api_authenticate(self):
        no_attempts, attempt, login_response = 5, 0, None
        while attempt < no_attempts:
            try:
                login_response = requests.post(
                    self.secrets.get_secret("IBF_API_URL") + "user/login",
                    data=[
                        ("email", self.secrets.get_secret("IBF_API_USER")),
                        ("password", self.secrets.get_secret("IBF_API_PASSWORD")),
                    ],
                )
                break
            except requests.exceptions.ConnectionError:
                attempt += 1
                logging.warning(
                    "IBF API currently not available, trying again in 1 minute"
                )
                time.sleep(60)
        if not login_response:
            raise ConnectionError("IBF API not available")
        return login_response.json()["user"]["token"]

    def ibf_api_post_request(self, path, body=None, files=None):
        token = self.__ibf_api_authenticate()
        if body is not None:
            headers = {
                "Authorization": "Bearer " + token,
                "Content-Type": "application/json",
                "Accept": "application/json",
            }
        elif files is not None:
            headers = {"Authorization": "Bearer " + token}
        else:
            raise ValueError("No body or files provided")
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        r = session.post(
            self.secrets.get_secret("IBF_API_URL") + path,
            json=body,
            files=files,
            headers=headers,
        )
        if r.status_code >= 400:
            raise ValueError(
                f"Error in IBF API POST request: {r.status_code}, {r.text}"
            )
        if not os.path.exists("logs"):
            os.makedirs("logs")
        if body:
            filename = body["date"]
            filename = "".join(x for x in filename if x.isalnum())
            filename = filename + ".json"
            filename = os.path.join("logs", filename)
            logs = {"endpoint": path, "payload": body}
            with open(filename, "a") as file:
                file.write(str(logs) + "\n")
        elif files:
            filename = datetime.today().strftime("%Y%m%d") + ".json"
            filename = os.path.join("logs", filename)
            logs = {"endpoint": path, "payload": files}
            with open(filename, "a") as file:
                file.write(str(logs) + "\n")

    def ibf_api_get_request(self, path, parameters=None):
        token = self.__ibf_api_authenticate()
        headers = {
            "Authorization": "Bearer " + token,
            "Accept": "*/*",
        }
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        r = session.get(
            self.secrets.get_secret("IBF_API_URL") + path,
            headers=headers,
            params=parameters,
        )
        if r.status_code >= 400:
            raise ValueError(f"Error in IBF API GET request: {r.status_code}, {r.text}")
        return r.json()

    def send_to_ibf_api(
        self,
        forecast_data: AdminDataSet,
        threshold_climateregion: ClimateRegionDataSet,
        forecast_climateregion: ClimateRegionDataSet,
        drought_extent: str = None,
        upload_time: datetime = datetime.now(),
    ):
        """Send drought forecast data to IBF API"""

        country = forecast_data.country
        # trigger_on_lead_time = self.settings.get_country_setting( country, "trigger-on-lead-time"    )
        admin_levels = self.settings.get_country_setting(country, "admin-levels")
        disasterType = "drought"
        indicators = [
            "population_affected",
            # "population_affected_percentage",
            "forecast_trigger",
            "forecast_severity",
        ]
        pipeline_will_trigger_portal = self.settings.get_country_setting(
            country, "pipeline-will-trigger-portal"
        )  # TODO: make varname more descriptive
        climate_regions_settings = self.settings.get_country_setting(
            country, "climate_region"
        )

        processed_pcodes = []

        current_year = upload_time.year
        current_month = upload_time.month
        current_month_abb = upload_time.strftime("%b")
        upload_time = datetime.today().strftime(
            f"{current_year}-{current_month:02}-%dT%H:%M:%SZ"
        )

        for climate_region_code in forecast_climateregion.get_climate_region_codes():
            climateregion = threshold_climateregion.get_data_unit(
                climate_region_code=climate_region_code
            )
            pcodes = climateregion.pcodes

            # get possible events for the climate region and current month
            possible_events = next(
                (
                    cr
                    for cr in climate_regions_settings
                    if cr["climate-region-code"] == climate_region_code
                )
            )["leadtime"][current_month_abb]

            lead_times_list = []
            expected_events = {}

            for entry in possible_events:
                for key, value in entry.items():
                    season_name = key
                    lead_time = int(value.split("-")[0])
                    expected_events[lead_time] = season_name
                    lead_times_list.append(lead_time)

            # current events
            events = self.__list_events_from_climateregion(
                forecast_climateregion, climate_region_code
            )

            climate_region_name = climateregion.climate_region_name

            for lead_time_event in range(0, 4):
                # NOTE: here we are assuming we will not expect two events in a climate region  with the same lead time
                if lead_time_event in list(expected_events.keys()):
                    season_name = expected_events[lead_time_event]
                    if climate_region_name.lower().split("_")[0] == "national":
                        event_name = f"{season_name}_National"
                    else:
                        event_name = (
                            f"{climate_region_name} {season_name}_{climate_region_name}"
                        )
                    # NOTE: exposure data is updated with new data during pre-season and not updated during the season
                    preseason_event, forecast_data_to_send = self.__fetch_or_fallback(
                        climate_region_code,
                        lead_time_event,
                        current_year,
                        current_month,
                        country,
                        season_name,
                        forecast_data,
                    )
                    if (preseason_event is False) or (
                        (lead_time_event not in events) and (preseason_event is None)
                    ):
                        # NOTE: this is to skip to upload empty exposure when no events and no data fetched
                        continue
                    for indicator in indicators:
                        for adm_level in admin_levels:
                            exposure_pcodes = []
                            for pcode in pcodes[f"{adm_level}"]:
                                forecast_admin = forecast_data_to_send.get_data_unit(
                                    pcode=pcode, lead_time=lead_time
                                )
                                amount = None
                                if indicator == "population_affected":
                                    amount = forecast_admin.pop_affected
                                elif indicator == "population_affected_percentage":
                                    amount = forecast_admin.pop_affected_perc
                                elif indicator == "forecast_severity":
                                    amount = forecast_admin.triggered
                                elif indicator == "forecast_trigger":
                                    amount = forecast_trigger_status(
                                        triggered=(forecast_admin.triggered >= 0),
                                        trigger_class=pipeline_will_trigger_portal,
                                    )
                                exposure_pcodes.append(
                                    {"placeCode": pcode, "amount": amount}
                                )
                                processed_pcodes.append(pcode)

                            body = {
                                "countryCodeISO3": country,
                                "leadTime": f"{lead_time_event}-month",
                                "dynamicIndicator": indicator,
                                "adminLevel": int(adm_level),
                                "exposurePlaceCodes": exposure_pcodes,
                                "disasterType": disasterType,
                                "eventName": event_name,
                                "date": upload_time,
                            }
                            statsPath = drought_extent.replace(
                                ".tif",
                                f"_{event_name}_{lead_time_event}-month_{country}_{adm_level}.json",
                            )
                            statsPath = statsPath.replace(
                                "rainfall_forecast", f"{indicator}"
                            )

                            with open(statsPath, "w") as fp:
                                json.dump(body, fp)

                            self.ibf_api_post_request(
                                "admin-area-dynamic-data/exposure", body=body
                            )
                    processed_pcodes = list(set(processed_pcodes))

        # END OF EVENT LOOP
        ###############################################################################################################

        # drought extent raster: admin-area-dynamic-data/raster/droughts
        self.rasters_sent = []

        for lead_time in range(0, 4):
            # NOTE: new drought extent raster is updated during the season
            drought_extent_new = drought_extent.replace(
                ".tif", f"_{lead_time}-month_{country}.tif"
            )

            # to accompdate file name requirement in IBF portal
            rainf_extent = drought_extent_new.replace(
                "rainfall_forecast", "rlower_tercile_probability"
            )
            rain_rp = drought_extent_new.replace("rainfall_forecast", "rain_rp")
            shutil.copy(
                rainf_extent, drought_extent_new.replace("rainfall_forecast", "rain_rp")
            )
            self.rasters_sent.append(rain_rp)
            files = {"file": open(rain_rp, "rb")}
            self.ibf_api_post_request(
                "admin-area-dynamic-data/raster/drought", files=files
            )

        # send empty exposure data
        if len(processed_pcodes) == 0:
            logging.info(f"send empty exposure data")
            for lead_time in set(lead_times_list):
                for indicator in indicators:
                    for adm_level in admin_levels:
                        exposure_pcodes = []
                        for pcode in forecast_data.get_pcodes(adm_level=adm_level):
                            if pcode not in processed_pcodes:
                                amount = None
                                if indicator == "population_affected":
                                    amount = 0
                                elif indicator == "population_affected_percentage":
                                    amount = 0.0
                                elif indicator == "forecast_trigger":
                                    amount = 0
                                elif indicator == "forecast_severity":
                                    amount = 0
                                exposure_pcodes.append(
                                    {"placeCode": pcode, "amount": amount}
                                )
                        body = {
                            "countryCodeISO3": country,
                            "leadTime": f"{lead_time}-month",  #  "1-day",  # this is a specific check IBF uses to establish no-trigger
                            "dynamicIndicator": indicator,
                            "adminLevel": adm_level,
                            "exposurePlaceCodes": exposure_pcodes,
                            "disasterType": disasterType,
                            "eventName": None,  # this is a specific check IBF uses to establish no-trigger
                            "date": upload_time,
                        }
                        self.ibf_api_post_request(
                            "admin-area-dynamic-data/exposure", body=body
                        )

                        statsPath = drought_extent.replace(
                            ".tif",
                            f"_null_{lead_time}-month_{country}_{adm_level}.json",
                        )
                        statsPath = statsPath.replace(
                            "rainfall_forecast", f"{indicator}"
                        )

                        with open(statsPath, "w") as fp:
                            json.dump(body, fp)

        # send notification
        body = {
            "countryCodeISO3": country,
            "disasterType": disasterType,
            "date": upload_time,
        }
        self.ibf_api_post_request("events/process", body=body)

    def save_pipeline_data(
        self, data_type: str, dataset: AdminDataSet, replace_country: bool = False
    ):
        """Upload pipeline datasets to Cosmos DB"""
        # To Cosmos DB
        if data_type not in COSMOS_DATA_TYPES:
            raise ValueError(
                f"Data type {data_type} is not supported."
                f"Supported storages are {', '.join(COSMOS_DATA_TYPES)}"
            )
        ## check data types
        if data_type == "seasonal-rainfall-forecast":
            for data_unit in dataset.data_units:
                if not isinstance(data_unit, ForecastDataUnit):
                    raise ValueError(
                        f"Data unit {data_unit} is not of type seasonal rainfall forecast"
                    )
        elif data_type == "seasonal-rainfall-forecast-climate-region":
            for data_unit in dataset.data_units:
                if not isinstance(data_unit, ForecastDataUnit):
                    raise ValueError(
                        f"Data unit {data_unit} is not of type seasonal rainfall forecast"
                    )
        elif data_type == "climate-region":
            for data_unit in dataset.data_units:
                if not isinstance(data_unit, ClimateRegionDataUnit):
                    raise ValueError(
                        f"Data unit {data_unit} is not of type ClimateregionDataUnit"
                    )

        client_ = cosmos_client.CosmosClient(
            self.secrets.get_secret("COSMOS_URL"),
            {"masterKey": self.secrets.get_secret("COSMOS_KEY")},
            user_agent="sml-api",
            user_agent_overwrite=True,
        )
        cosmos_db = client_.get_database_client("drought-pipeline")
        cosmos_container_client = cosmos_db.get_container_client(data_type)

        # Check if data of the same month exists
        start_date, end_date = self.__dates_for_month(
            current_year=dataset.timestamp.year,
            current_month=dataset.timestamp.month,
        )
        if replace_country:
            query = get_cosmos_query(
                start_date=start_date, end_date=end_date, country=dataset.country
            )
            old_records = cosmos_container_client.query_items(query)
            for old_record in old_records:
                cosmos_container_client.delete_item(
                    item=old_record.get("id"), partition_key=dataset.country
                )

        # Upsert new records
        for data_unit in dataset.data_units:
            record = vars(data_unit)
            record["timestamp"] = dataset.timestamp.strftime("%Y-%m-%dT%H:%M:%S")
            record["country"] = dataset.country
            record["id"] = get_data_unit_id(data_unit, dataset)
            cosmos_container_client.upsert_item(body=record)

    def get_pipeline_data(
        self,
        data_type,
        country,
        start_date=None,
        end_date=None,
        adm_level=None,
        pcode=None,
        lead_time=None,
    ) -> AdminDataSet:
        """Download pipeline datasets from Cosmos DB"""
        if data_type not in COSMOS_DATA_TYPES:
            raise ValueError(
                f"Data type {data_type} is not supported."
                f"Supported storages are {', '.join(COSMOS_DATA_TYPES)}"
            )
        client_ = cosmos_client.CosmosClient(
            self.secrets.get_secret("COSMOS_URL"),
            {"masterKey": self.secrets.get_secret("COSMOS_KEY")},
            user_agent="ibf-flood-pipeline",
            user_agent_overwrite=True,
        )
        cosmos_db = client_.get_database_client("drought-pipeline")
        cosmos_container_client = cosmos_db.get_container_client(data_type)
        query = get_cosmos_query(
            start_date, end_date, country, adm_level, pcode, lead_time
        )
        records_query = cosmos_container_client.query_items(
            query=query,
            enable_cross_partition_query=(
                True if country is None else None
            ),  # country must be the partition key
        )
        records = []
        for record in records_query:
            records.append(copy.deepcopy(record))
        datasets = []
        countries = list(set([record["country"] for record in records]))
        timestamps = list(set([record["timestamp"] for record in records]))
        for country in countries:
            for timestamp in timestamps:
                data_units = []
                for record in records:
                    if (
                        record["country"] == country
                        and record["timestamp"] == timestamp
                    ):
                        if data_type in [
                            "seasonal-rainfall-forecast",
                            "seasonal-rainfall-forecast-climate-region",
                        ]:
                            data_unit = ForecastDataUnit(
                                adm_level=record["adm_level"],
                                pcode=record["pcode"],
                                climate_region_code=record["climate_region_code"],
                                lead_time=record["lead_time"],
                                triggered=record["triggered"],
                                tercile_upper=record["tercile_upper"],
                                tercile_lower=record["tercile_lower"],
                                likelihood=record["likelihood"],
                                pop_affected=record["pop_affected"],
                                pop_affected_perc=record["pop_affected_perc"],
                                alert_class=record["alert_class"],
                            )
                        # elif data_type == "seasonal-rainfall-forecast-climate-region": # TODO:refine
                        #     data_unit = ForecastDataUnit(
                        #         adm_level=record["adm_level"],
                        #         pcode=record["pcode"],
                        #         lead_time=record["lead_time"],
                        #         triggered=record["triggered"],
                        #         tercile_upper=record["tercile_upper"],
                        #         tercile_lower=record["tercile_lower"],
                        #         likelihood=record["likelihood"],
                        #         pop_affected=record["pop_affected"],
                        #         pop_affected_perc=record["pop_affected_perc"],
                        #         alert_class=record["alert_class"],
                        #     )
                        elif data_type == "climate-region":
                            data_unit = ClimateRegionDataUnit(
                                adm_level=record["adm_level"],
                                climate_region_code=record["climate_region_code"],
                                climate_region_name=record["climate_region_name"],
                                pcodes=record["pcodes"],
                            )

                        else:
                            raise ValueError(f"Invalid data type {data_type}")
                        data_units.append(data_unit)
                if data_type in [
                    "seasonal-rainfall-forecast"
                ]:  # , "seasonal-rainfall-forecast-climate-region"]
                    adm_levels = list(
                        set([data_unit.adm_level for data_unit in data_units])
                    )
                    dataset = AdminDataSet(
                        country=country,
                        timestamp=timestamp,
                        adm_levels=adm_levels,
                        data_units=data_units,
                    )
                    datasets.append(dataset)
                else:
                    dataset = ClimateRegionDataSet(
                        country=country,
                        timestamp=timestamp,
                        data_units=data_units,
                    )
                    datasets.append(dataset)
        if len(datasets) == 0:
            raise KeyError(
                f"No datasets of type '{data_type}' found for country {country} in date range "
                f"{start_date} - {end_date}."
            )
        elif len(datasets) > 1:
            logging.warning(
                f"Multiple datasets of type '{data_type}' found for country {country} in date range "
                f"{start_date} - {end_date}; returning the latest (timestamp {datasets[-1].timestamp}). "
            )
        return datasets[-1]

    def __get_blob_service_client(self, blob_path: str):
        """Get service client for Azure Blob Storage"""
        blob_service_client = BlobServiceClient.from_connection_string(
            f"DefaultEndpointsProtocol=https;"
            f'AccountName={self.secrets.get_secret("BLOB_ACCOUNT_NAME")};'
            f'AccountKey={self.secrets.get_secret("BLOB_ACCOUNT_KEY")};'
            f"EndpointSuffix=core.windows.net"
        )
        container = self.settings.get_setting("blob_container")
        # blob_path = self.settings.get_setting("blob_storage_path")
        return blob_service_client.get_blob_client(container=container, blob=blob_path)

    def save_to_blob(self, local_path: str, file_dir_blob: str):
        """Save file to Azure Blob Storage"""
        # upload to Azure Blob Storage
        blob_client = self.__get_blob_service_client(file_dir_blob)
        with open(local_path, "rb") as upload_file:
            blob_client.upload_blob(upload_file, overwrite=True)

    def upload_json_files(self, country, local_path: str):
        """Find all JSON files in a directory and upload them to Azure Blob Storage"""
        if not os.path.exists(local_path):
            logging.info(f"Directory not found: {local_path}")
            return

        # Create a new folder in Blob Storage based on today's date
        blob_storage_path = self.settings.get_setting("databases")["blob_storage_path"]
        month = datetime.now().strftime("%Y-%m")
        blob_folder = os.path.join(
            blob_storage_path, country, month
        )  # if blob_storage_path else today_date

        for file_name in os.listdir(local_path):
            if file_name.endswith(".json") or file_name.endswith(
                ".tif"
            ):  #            if file_name.endswith(".json"):
                local_path_ = os.path.join(local_path, file_name)
                blob_path = os.path.join(
                    blob_folder, file_name
                )  # if blob_folder else file_name
                logging.info(
                    f"Uploading {local_path_} to {blob_path} in Blob Storage..."
                )
                self.save_to_blob(local_path_, blob_path)

    def get_from_blob(self, local_path: str, blob_path: str):
        """Get file from Azure Blob Storage"""
        blob_client = self.__get_blob_service_client(blob_path)

        with open(local_path, "wb") as download_file:
            try:
                download_file.write(blob_client.download_blob().readall())
            except ResourceNotFoundError:
                raise FileNotFoundError(
                    f"File {blob_path} not found in Azure Blob Storage"
                )

    def download_ecmwf_forecast(self, country, data_dir, current_year, current_month):
        """Download ECMWF seasonal hindcast data for historical period
        Args:
            country (str): Country name
            data_dir (str): Directory to save data
            current_year (int): Current year
            current_month (int): Current month
        """
        gdf = self.get_adm_boundaries(country, 1)

        min_x, min_y, max_x, max_y = gdf.total_bounds

        KEY = os.getenv("CDSAPI_KEY")
        URL = "https://cds.climate.copernicus.eu/api"

        c = cdsapi.Client(url=URL, key=KEY, wait_until_complete=False, delete=False)

        # Forecast data request
        dataset = "seasonal-monthly-single-levels"
        request = {
            "originating_centre": "ecmwf",
            "system": "51",
            "variable": ["total_precipitation"],
            "product_type": ["monthly_mean"],
            "year": [current_year],
            "month": [current_month],
            "leadtime_month": ["1", "2", "3", "4", "5", "6"],
            "data_format": "grib",
            "area": [
                int(x) for x in [max_y + 1, min_x - 1, min_y - 1, max_x + 1]
            ],  # North, West, South, East
        }
        target = f"{data_dir}/ecmwf_seas5_forecast_monthly_tp.grib"
        c.retrieve(dataset, request, target)

        sleep = 30
        time.sleep(sleep)

        request = {
            "originating_centre": "ecmwf",
            "system": "51",
            "variable": ["total_precipitation"],
            "product_type": ["monthly_mean"],
            "year": [
                "1991",
                "1992",
                "1993",
                "1994",
                "1995",
                "1996",
                "1997",
                "1998",
                "1999",
                "2000",
                "2001",
                "2002",
                "2003",
                "2004",
                "2005",
                "2006",
                "2007",
                "2008",
                "2009",
                "2010",
                "2011",
                "2012",
                "2013",
                "2014",
                "2015",
                "2016",
                "2017",
                "2018",
                "2019",
                "2020",
            ],
            "month": ["03"],
            "leadtime_month": ["1", "2", "3", "4", "5", "6"],
            "data_format": "grib",
            "area": [
                int(x) for x in [max_y + 1, min_x - 1, min_y - 1, max_x + 1]
            ],  # North, West, South, East
        }
        target = f"{data_dir}/ecmwf_seas5_hindcast_monthly_tp.grib"
        c.retrieve(dataset, request, target)

    def __look_up_dates(
        self,
        country,
        climate_region_code,
        season_name,
        current_year,
        current_month,
        lead_time_event="1-month",
    ):
        """Look up the month of the most recent 1-month forecast
        and return the month's start and end date.
        Args:
            country (str): Country name
            event_name (str): Event name
            lead_time_event (int): Lead time event"""
        datestart, dateend = None, None
        climate_regions_settings = self.settings.get_country_setting(
            country, "climate_region"
        )
        # get the leadtime month dictionary for the climate region
        month_dict = next(
            (
                cr
                for cr in climate_regions_settings
                if cr["climate-region-code"] == climate_region_code
            )
        )["leadtime"]
        # Iterate in reverse to get the most recent one
        for month in reversed(list(month_dict.keys())):
            for entry in month_dict[month]:
                if season_name in entry and entry[season_name] == lead_time_event:
                    month_num = datetime.strptime(month, "%b").month
                    if month_num > current_month:
                        current_year -= 1
                    datestart, dateend = self.__dates_for_month(
                        current_year=current_year,
                        current_month=month_num,
                    )
                    return datestart, dateend
        if datestart is None or dateend is None:
            raise ValueError("No matching forecast period found.")
        return datestart.date(), dateend.date()

    def __dates_for_month(self, current_year, current_month):
        """Return start and end date of the given year and month."""
        datestart = datetime(current_year, current_month, 1)
        if current_month == 12:
            next_month = datetime(current_year + 1, 1, 1)
        else:
            next_month = datetime(current_year, current_month + 1, 1)
        dateend = next_month - timedelta(days=1)
        return datestart.date(), dateend.date()

    def __fetch_or_fallback(
        self,
        climate_region_code,
        lead_time,
        current_year,
        current_month,
        country,
        season_name,
        forecast_data,
    ):
        """
        Fetch data of the right events and lead time from the database or return fallback data if not available.
        """
        if lead_time == 0:
            datestart, dateend = self.__look_up_dates(
                country,
                climate_region_code,
                season_name,
                current_year=current_year,
                current_month=current_month,
                lead_time_event="1-month",
            )
            try:
                logging.info(
                    f"fetch climate region {climate_region_code} between datestart: {datestart}, dateend: {dateend}"
                )
                preseason_forecast_climateregion = self.get_pipeline_data(
                    "seasonal-rainfall-forecast-climate-region",
                    country=country,
                    start_date=datestart,
                    end_date=dateend,
                )
                preseason_events = self.__list_events_from_climateregion(
                    preseason_forecast_climateregion, climate_region_code
                )
                if lead_time in list(preseason_events.keys()):
                    preseason_event_lead_time_0 = (
                        True if preseason_events[0] == "trigger" else False
                    )
                else:
                    preseason_event_lead_time_0 = False
                logging.info(
                    f"fetch forecast data between datestart: {datestart}, dateend: {dateend}"
                )
                forecast_data = self.get_pipeline_data(
                    "seasonal-rainfall-forecast",
                    country=country,
                    start_date=datestart,
                    end_date=dateend,
                )
            except Exception as e:
                logging.warning(
                    f"Fetching data failed. {e}. Fallback to current forecast data."
                )
                preseason_event_lead_time_0 = None
        else:
            preseason_event_lead_time_0 = None
        return preseason_event_lead_time_0, forecast_data

    def __list_events_from_climateregion(
        self, forecast_climateregion, climate_region_code
    ):
        """
        List foreacast events from climateregion data unit
        """
        events = {}
        # triggered_lead_times = []
        for lead_time in range(0, 6):
            if (
                forecast_climateregion.get_climate_region_data_unit(
                    climate_region_code, lead_time
                ).alert_class
                != "no"
            ):
                events[lead_time] = "alert"
        for lead_time in range(0, 6):
            if forecast_climateregion.get_climate_region_data_unit(
                climate_region_code, lead_time
            ).triggered:
                events[lead_time] = "trigger"
                # triggered_lead_times.append(lead_time)
        events = dict(sorted(events.items()))
        return events
