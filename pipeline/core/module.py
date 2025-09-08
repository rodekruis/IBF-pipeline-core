from pipeline.core.load import Load
from pipeline.core.secrets import Secrets
from pipeline.core.settings import Settings
from pipeline.core.data import PipelineDataSets


class PipelineModule:
    """Base class for pipeline module"""

    def __init__(self, settings: Settings, secrets: Secrets, country: str):
        self.settings = settings
        if country not in [c["name"] for c in self.settings.get_setting("countries")]:
            raise ValueError(f"No config found for country {country}")
        self.country = country
        self.data = PipelineDataSets(country=country)
        self.load = Load(country=country, settings=settings, secrets=secrets)
