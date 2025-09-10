from ibfpipeline.riverflood.extract import Extract
from ibfpipeline.riverflood.forecast import Forecast
from ibfpipeline.riverflood.load import RiverFloodLoad
from ibfpipeline.riverflood.data import RiverFloodPipelineDataSets
from ibfpipeline.core.secrets import Secrets
from ibfpipeline.core.settings import Settings
from ibfpipeline.core.logger import logger


class Pipeline:
    """Base class for flood data pipeline"""

    def __init__(self, settings: Settings, secrets: Secrets, country: str):
        self.settings = settings
        if country not in [c["name"] for c in self.settings.get_setting("countries")]:
            raise ValueError(f"No config found for country {country}")
        self.country = country

        # Initialize empty data sets
        self.data = RiverFloodPipelineDataSets(country=country, settings=settings)

        # Initialize pipeline modules
        self.load = RiverFloodLoad(
            country=country, settings=settings, secrets=secrets, data=self.data
        )
        self.extract = Extract(
            country=country,
            settings=settings,
            secrets=secrets,
            data=self.data,
        )
        self.forecast = Forecast(
            country=country,
            settings=settings,
            secrets=secrets,
            data=self.data,
        )

        # Load thresholds
        self.data.threshold_admin = self.load.get_thresholds_admin()
        self.data.threshold_station = self.load.get_thresholds_station()

    def run_pipeline(
        self,
        prepare: bool = True,
        forecast: bool = True,
        send: bool = True,
        debug: bool = False,  # fast extraction on yesterday's data, using only one ensemble member
    ):
        """Run the flood data pipeline"""

        if prepare:
            logger.info("prepare discharge data")
            self.extract.prepare_glofas_data(country=self.country, debug=debug)

        if forecast:
            logger.info(f"extract discharge data")
            self.extract.extract_glofas_data(country=self.country, debug=debug)
            logger.info("forecast floods")
            self.forecast.compute_forecast()

        if send:
            logger.info("send data to IBF API")
            self.load.send_to_ibf_api(
                flood_extent=self.forecast.flood_extent_raster,
            )
