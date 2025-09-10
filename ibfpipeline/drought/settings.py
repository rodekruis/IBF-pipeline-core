import os
import yaml
from urllib.parse import urlparse
from typing import List
from datetime import date

# Default config file path
DEFAULT_CONFIG_PATH = os.path.join(os.path.dirname(__file__), "./config/config.yml")


def is_url(x):
    try:
        result = urlparse(x)
        return all([result.scheme, result.netloc])
    except ValueError:
        return False


class Settings:
    """
    Settings (URLs, constants, etc.)
    """

    def __init__(self, path_or_url: str):  # = DEFAULT_CONFIG_PATH):
        self.setting_source: str = "yaml"
        self.setting_path: str = path_or_url
        if not os.path.exists(self.setting_path) and not is_url(self.setting_path):
            raise ValueError(f"Setting file {self.setting_path} not found")
        with open(self.setting_path, "r") as file:
            self.settings = yaml.load(file, Loader=yaml.FullLoader)

    def get_setting(self, setting: str):
        setting_value = None
        if setting in self.settings.keys():
            setting_value = self.settings[setting]
        else:
            for key in self.settings.keys():
                if type(self.settings[key]) == dict:
                    if setting in self.settings[key].keys():
                        setting_value = self.settings[key][setting]
                elif type(self.settings[key]) == list:
                    for i in range(len(self.settings[key])):
                        if setting in self.settings[key][i].keys():
                            setting_value = self.settings[key][i][setting]
        if not setting_value:
            raise ValueError(f"Setting {setting} not found in {self.setting_path}")
        return setting_value

    def get_country_setting(self, country: str, setting: str):
        country_setting = next(
            x for x in self.get_setting("countries") if x["name"] == country
        )
        if not country_setting:
            raise ValueError(f"Country {country} not found in {self.setting_path}")
        if setting not in country_setting.keys():
            raise ValueError(f"Setting {setting} not found for country {country}")
        return country_setting[setting]

    def check_settings(self, settings: List[str]):
        missing_settings = []
        for setting in settings:
            try:
                self.get_setting(setting)
            except ValueError:
                missing_settings.append(setting)
        if missing_settings:
            raise Exception(
                f"Missing settings {', '.join(missing_settings)} in {self.setting_path}"
            )

    # def get_leadtime_for_climate_region(self, country: str, climate_region: str, month: str):
    #     # Fetch country setting
    #     country_setting = next(
    #         (x for x in self.get_setting("countries") if x["name"] == country), None
    #     )
    #
    #     if not country_setting:
    #         raise ValueError(f"Country {country} not found in {self.setting_path}")
    #
    #     # Loop through the climate regions for the country
    #     for region in country_setting.get("Climate_Region", []):
    #         if region["name"] == climate_region:
    #             # Check if the month exists in the leadtime for this climate region
    #             leadtime = region.get("leadtime", {}).get(month)
    #             if leadtime:
    #                 return leadtime
    #             else:
    #                 raise ValueError(f"Month {month} not found for climate region {climate_region} in country {country}")
    #
    #     raise ValueError(f"Climate region {climate_region} not found in country {country}")
    #
    # def get_leadtime_for_climate_region_code(self, country: str, climate_region_code: int, month: str):
    #     # Fetch country setting
    #     country_setting = next(
    #         (x for x in self.get_setting("countries") if x["name"] == country), None
    #     )
    #
    #     if not country_setting:
    #         raise ValueError(f"Country {country} not found in {self.setting_path}")
    #
    #     # Loop through the climate regions for the country
    #     for region in country_setting.get("climate_region", []):
    #         if region["climate-region-code"] == climate_region_code:
    #             # Check if the month exists in the leadtime for this climate region
    #             leadtime = region.get("leadtime", {}).get(month)
    #             if leadtime:
    #                 return leadtime
    #             else:
    #                 raise ValueError(f"Month {month} not found for climate region code {climate_region_code} in country {country}")
    #
    #     raise ValueError(f"Climate region code {climate_region_code} not found in country {country}")
    #
    # def get_climate_region_name_by_code(self, country: str, climate_region_code: int):
    #     # Fetch country setting
    #     country_setting = next(
    #         (x for x in self.get_setting("countries") if x["name"] == country), None
    #     )
    #
    #     if not country_setting:
    #         raise ValueError(f"Country {country} not found in {self.setting_path}")
    #
    #     # Loop through the climate regions for the country
    #     for region in country_setting.get("climate_region", []):
    #         if region["climate-region-code"] == climate_region_code:
    #             # Return the name of the climate region for the matched code
    #             return region["name"]
    #
    #     raise ValueError(f"Climate region code {climate_region_code} not found in country {country}")
    #
    #
    # def get_all_leadtime_for_climate_region_code(self, country: str, climate_region_code: int):
    #     # Fetch country setting
    #     country_setting = next(
    #         (x for x in self.get_setting("countries") if x["name"] == country), None
    #     )
    #
    #     if not country_setting:
    #         raise ValueError(f"Country {country} not found in {self.setting_path}")
    #
    #     # Loop through the climate regions for the country
    #     for region in country_setting.get('climate_region', []):
    #         if region["climate-region-code"] == climate_region_code:
    #             return region['leadtime']
    #
    #     raise ValueError(f"Climate region code {climate_region_code} not found in country {country}")
