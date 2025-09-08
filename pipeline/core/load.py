from __future__ import annotations

import time

from pipeline.core.secrets import Secrets
from pipeline.core.settings import Settings
from pipeline.core.logger import logger

from urllib.error import HTTPError
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import requests
import geopandas as gpd
from azure.storage.blob import BlobServiceClient
from azure.core.exceptions import ResourceNotFoundError


class Load:
    """Download/upload data from/to a data storage"""

    def __init__(self, country: str, settings: Settings, secrets: Secrets):
        self.country = country
        self.settings = self.check_settings(settings)
        self.secrets = self.check_secrets(secrets)
        self.rasters_sent = []

    def check_settings(self, settings: Settings):
        """Check settings"""
        if not isinstance(settings, Settings):
            raise TypeError(f"invalid format of settings, use settings.Settings")
        settings.check_settings(
            ["postgresql_server", "postgresql_port", "postgresql_database"]
        )
        return settings

    def check_secrets(self, secrets: Secrets):
        """Check secrets for storage"""
        if not isinstance(secrets, Secrets):
            raise TypeError(f"invalid format of secrets, use secrets.Secrets")
        secrets.check_secrets(
            [
                "BLOB_ACCOUNT_NAME",
                "BLOB_ACCOUNT_KEY",
                "IBF_API_URL",
                "IBF_API_USER",
                "IBF_API_PASSWORD",
            ]
        )
        return secrets

    def get_population_density(self, file_path: str):
        """Get population density data from worldpop and save to file_path"""
        r = requests.get(
            f"{self.settings.get_setting('worldpop_url')}/{self.country}/{self.country.lower()}_ppp_2022_1km_UNadj_constrained.tif"
        )
        if "404 Not Found" in str(r.content):
            raise FileNotFoundError(
                f"Population density data not found for country {self.country}"
            )
        with open(file_path, "wb") as file:
            file.write(r.content)

    def get_adm_boundaries(self, adm_level: int) -> gpd.GeoDataFrame:
        """Get admin areas from IBF API"""
        try:
            adm_boundaries = self.ibf_api_request(
                "GET",
                f"admin-areas/{self.country}/{adm_level}",
            )
        except HTTPError:
            raise FileNotFoundError(
                f"Admin areas for country {self.country}"
                f" and admin level {adm_level} not found"
            )
        gdf_adm_boundaries = gpd.GeoDataFrame.from_features(adm_boundaries["features"])
        gdf_adm_boundaries.set_crs(epsg=4326, inplace=True)
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
                logger.warning(
                    "IBF API currently not available, trying again in 1 minute"
                )
                time.sleep(60)
        if not login_response:
            raise ConnectionError("IBF API not available")
        return login_response.json()["user"]["token"]

    def ibf_api_request(self, method, path, parameters=None, body=None, files=None):
        token = self.__ibf_api_authenticate()
        headers = {"Authorization": "Bearer " + token}
        if parameters is not None:
            headers["Accept"] = "*/*"
        elif body is not None:
            headers["Content-Type"] = "application/json"
            headers["Accept"] = "application/json"
        session = requests.Session()
        retry = Retry(connect=3, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        if method.upper() == "POST":
            r = session.post(
                self.secrets.get_secret("IBF_API_URL") + path,
                json=body,
                files=files,
                headers=headers,
            )
        elif method.upper() == "GET":
            r = session.get(
                self.secrets.get_secret("IBF_API_URL") + path,
                headers=headers,
                params=parameters,
            )
        else:
            raise KeyError(f"Method {method} is not supported for IBF API")
        if r.status_code >= 400:
            raise ValueError(
                f"Error in IBF API POST request: {r.status_code}, {r.text}"
            )

    def __get_blob_service_client(self, blob_path: str):
        """Get service client for Azure Blob Storage"""
        blob_service_client = BlobServiceClient.from_connection_string(
            f"DefaultEndpointsProtocol=https;"
            f'AccountName={self.secrets.get_secret("BLOB_ACCOUNT_NAME")};'
            f'AccountKey={self.secrets.get_secret("BLOB_ACCOUNT_KEY")};'
            f"EndpointSuffix=core.windows.net"
        )
        container = self.settings.get_setting("blob_container")
        return blob_service_client.get_blob_client(container=container, blob=blob_path)

    def save_to_blob(self, local_path: str, file_dir_blob: str):
        """Save file to Azure Blob Storage"""
        # upload to Azure Blob Storage
        blob_client = self.__get_blob_service_client(file_dir_blob)
        with open(local_path, "rb") as upload_file:
            blob_client.upload_blob(upload_file, overwrite=True)

    def get_from_blob(self, local_path: str, blob_path: str):
        """Get file from Azure Blob Storage"""
        blob_client = self.__get_blob_service_client(blob_path)

        with open(local_path, "wb") as download_file:
            try:
                download_file.write(blob_client.download_blob(timeout=120).readall())
            except ResourceNotFoundError:
                raise FileNotFoundError(
                    f"File {blob_path} not found in Azure Blob Storage"
                )
