from datetime import datetime
from typing import List
from ibfpipeline.core.settings import Settings
import numpy as np


class AdminDataUnit:
    """Base class for admin data units"""

    def __init__(self, **kwargs):
        self.adm_level: int = kwargs.get("adm_level")
        self.pcode: str = kwargs.get("pcode")


class ClimateRegionDataUnit:
    """Base class for climate region data units"""

    def __init__(self, **kwargs):
        self.climate_region_code: str = kwargs.get("climate_region_code")
        self.climate_region_name: str = kwargs.get("climate_region_name")
        self.adm_level: int = kwargs.get("adm_level")
        self.pcodes: dict = kwargs.get(
            "pcodes"
        )  # pcodes of associated administrative divisions


class RainfallDataUnit(AdminDataUnit):
    """Rainfall data unit - admin"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.lead_time: int = kwargs.get("lead_time", 0)
        self.tercile_lower: float = kwargs.get("tercile_lower", None)
        self.tercile_upper: float = kwargs.get("tercile_upper", None)
        self.rainfall_forecast: List[float] = kwargs.get("rainfall_forecast", None)
        self.likelihood: float = kwargs.get("likelihood", None)
        self.trigger: bool = kwargs.get("trigger", None)

        if hasattr(self.likelihood, "__iter__"):
            self.compute_threshold()

    def compute_threshold(self):
        """Compute the percentage of forecast values below the tercile lower threshold"""
        self.likelihood = {}
        for lead_time in range(1, 7):
            # Check if necessary data is available
            if lead_time not in self.rainfall_forecast:
                print(
                    f"Warning: Missing rainfall forecast data for lead_time {lead_time}"
                )  # TODO: replace all print() with logger
                continue
            if lead_time not in self.tercile_lower:
                print(f"Warning: Missing tercile lower data for lead_time {lead_time}")
                continue

            forecast_values = np.array(self.rainfall_forecast[lead_time])
            tercile_lower_value = self.tercile_lower[lead_time]

            # Ensure forecast values are not empty
            if len(forecast_values) == 0:
                print(
                    f"Warning: Empty rainfall forecast values for lead_time {lead_time}"
                )
                continue

            # Compute the percentage of values below the tercile lower threshold
            percentage_below_tercile_lower = (
                forecast_values < tercile_lower_value
            ).sum() / len(forecast_values)

            key = f"{lead_time}_month"
            self.likelihood[key] = int(percentage_below_tercile_lower)


class RainfallClimateRegionDataUnit(ClimateRegionDataUnit):
    """rainfall data unit - climate region"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.lead_time: int = kwargs.get("lead_time", 0)
        self.tercile_lower: float = kwargs.get("tercile_lower", None)
        self.tercile_upper: float = kwargs.get("tercile_upper", None)
        self.rainfall_forecast: List[float] = kwargs.get("rainfall_forecast", None)
        self.likelihood: float = kwargs.get("likelihood", None)
        self.trigger: bool = kwargs.get("trigger", None)

        if hasattr(self.likelihood, "__iter__"):
            self.compute_threshold()

    def compute_threshold(self):
        """Compute the percentage of forecast values below the tercile lower threshold"""
        self.likelihood = {}
        for lead_time in range(1, 7):
            # Check if necessary data is available
            if lead_time not in self.rainfall_forecast:
                print(
                    f"Warning: Missing rainfall forecast data for lead_time {lead_time}"
                )
                continue
            if lead_time not in self.tercile_lower:
                print(f"Warning: Missing tercile lower data for lead_time {lead_time}")
                continue

            forecast_values = np.array(self.rainfall_forecast[lead_time])
            tercile_lower_value = self.tercile_lower[lead_time]

            # Ensure forecast values are not empty
            if len(forecast_values) == 0:
                print(
                    f"Warning: Empty rainfall forecast values for lead_time {lead_time}"
                )
                continue

            # Compute the percentage of values below the tercile lower threshold
            percentage_below_tercile_lower = (
                forecast_values < tercile_lower_value
            ).sum() / len(forecast_values)

            key = f"{lead_time}_month"
            self.likelihood[key] = int(percentage_below_tercile_lower)


class ForecastDataUnit(AdminDataUnit):
    """Drought forecast data unit"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.lead_time: int = kwargs.get("lead_time")

        # Expecting forecasts as a dictionary with 'tercile_lower', 'tercile_upper', and 'forecast' keys
        self.tercile_lower: float = kwargs.get("tercile_lower", None)
        self.tercile_upper: float = kwargs.get("tercile_upper", None)
        self.forecast: list = kwargs.get("forecast", None)
        self.season: str = kwargs.get("season", None)
        self.climate_region_code: int = kwargs.get("climate_region_code", None)
        self.climate_region_name: str = kwargs.get("climate_region_name", None)
        self.pop_affected: int = kwargs.get("pop_affected", 0)
        self.pop_affected_perc: float = kwargs.get("pop_affected_perc", 0.0)
        self.triggered: bool = kwargs.get("triggered", None)
        self.likelihood: float = kwargs.get("likelihood", None)
        self.return_period: float = kwargs.get("return_period", None)
        self.alert_class: str = kwargs.get("alert_class", None)


class AdminDataSet:
    """Base class for admin data sets"""

    def __init__(
        self,
        country: str = None,
        timestamp: datetime = datetime.now(),
        adm_levels: List[int] = None,
        data_units: List[AdminDataUnit] = None,
    ):
        self.country = country
        self.timestamp = timestamp
        self.adm_levels = adm_levels
        self.data_units = data_units

    def get_pcodes(self, adm_level: int = None):
        """Return list of unique pcodes, optionally filtered by adm_level"""
        if not adm_level:
            return list(set([x.pcode for x in self.data_units]))
        else:
            return list(
                set([x.pcode for x in self.data_units if x.adm_level == adm_level])
            )

    def get_climate_region_codes(self):
        """Return list of unique climate_region_code"""
        return list(
            set(
                [
                    x.climate_region_code
                    for x in self.data_units
                    if hasattr(x, "climate_region_code")
                ]
            )
        )

    def get_lead_times(self):
        """Return list of unique lead times"""
        return list(
            set([x.lead_time for x in self.data_units if hasattr(x, "lead_time")])
        )

    def get_data_units(self, lead_time: int = None, adm_level: int = None):
        """Return list of data units filtered by lead time and/or admin level"""
        if not self.data_units:
            raise ValueError("Data units not found")
        if lead_time is not None and adm_level is not None:
            return list(
                filter(
                    lambda x: x.lead_time == lead_time and x.adm_level == adm_level,
                    self.data_units,
                )
            )
        elif lead_time is not None:
            return list(filter(lambda x: x.lead_time == lead_time, self.data_units))
        elif adm_level is not None:
            return list(filter(lambda x: x.adm_level == adm_level, self.data_units))
        else:
            return self.data_units

    def get_data_unit(self, pcode: str, lead_time: int = None) -> AdminDataUnit:
        """Get data unit by pcode and optionally by lead time"""
        if not self.data_units:
            raise ValueError("Data units not found")
        if lead_time is not None:
            bdu = next(
                filter(
                    lambda x: x.pcode == pcode and x.lead_time == lead_time,
                    self.data_units,
                ),
                None,
            )
        else:
            bdu = next(
                filter(lambda x: x.pcode == pcode, self.data_units),
                None,
            )
        if not bdu:
            raise ValueError(
                f"Data unit with pcode {pcode} and lead_time {lead_time} not found"
            )
        else:
            return bdu

    def get_data_unit_climate_region(
        self, climate_region_code: str, lead_time: int = None
    ) -> AdminDataUnit:
        """Get data unit by pcode and optionally by lead time"""
        if not self.data_units:
            raise ValueError("Data units not found")
        if lead_time is not None:
            bdu = next(
                filter(
                    lambda x: x.climate_region_code == climate_region_code
                    and x.lead_time == lead_time,
                    self.data_units,
                ),
                None,
            )
        else:
            bdu = next(
                filter(
                    lambda x: x.climate_region_code == climate_region_code,
                    self.data_units,
                ),
                None,
            )
        if not bdu:
            raise ValueError(
                f"Data unit with pcode {climate_region_code} and lead_time {lead_time} not found"
            )
        else:
            return bdu

    def upsert_data_unit(self, data_unit: AdminDataUnit):
        """Add data unit; if it already exists, update it"""
        if not self.data_units:
            self.data_units = [data_unit]
        if hasattr(data_unit, "lead_time"):
            bdu = next(
                filter(
                    lambda x: x[1].pcode == data_unit.pcode
                    and x[1].lead_time == data_unit.lead_time,
                    enumerate(self.data_units),
                ),
                None,
            )
        else:
            bdu = next(
                filter(
                    lambda x: x[1].pcode == data_unit.pcode,
                    enumerate(self.data_units),
                ),
                None,
            )
        if not bdu:
            self.data_units.append(data_unit)
        else:
            self.data_units[bdu[0]] = data_unit

    def is_any_triggered(self):
        """Check if any data unit is triggered"""
        if not self.data_units:
            raise ValueError("Data units not found")
        if type(self.data_units[0]) != ForecastDataUnit:
            raise ValueError("Data units are not forecast data units")
        return any([x.triggered for x in self.data_units])


class ClimateRegionDataSet:
    """Base class for station data sets"""

    def __init__(
        self,
        country: str = None,
        timestamp: datetime = datetime.now(),
        data_units: List["ClimateRegionDataUnit"] = None,
    ):
        self.country = country
        self.timestamp = timestamp
        self.data_units = data_units

    def get_data_unit(self, climate_region_code: str) -> "ClimateRegionDataUnit":
        """Get data unit by climate_region_code"""
        if not self.data_units:
            raise ValueError("Data units not found")

        bdu = next(
            filter(
                lambda x: x.climate_region_code == climate_region_code, self.data_units
            ),
            None,
        )

        if not bdu:
            raise ValueError(
                f"Data unit with climate_region_code {climate_region_code} not found"
            )
        return bdu

    def get_climate_region_data_unit(
        self, climate_region_code: str, lead_time: int = None
    ) -> AdminDataUnit:
        """Get data unit by climate_region_code and optionally by lead time"""
        if not self.data_units:
            raise ValueError("Data units not found")
        if lead_time is not None:
            bdu = next(
                filter(
                    lambda x: x.climate_region_code == climate_region_code
                    and x.lead_time == lead_time,
                    self.data_units,
                ),
                None,
            )
        else:
            bdu = next(
                filter(
                    lambda x: x.climate_region_code == climate_region_code,
                    self.data_units,
                ),
                None,
            )
        if not bdu:
            raise ValueError(
                f"Data unit with climate_region_code {climate_region_code} and lead_time {lead_time} not found"
            )
        else:
            return bdu

    def get_data_units(self, lead_time: int = None, adm_level: int = None):
        """Return list of data units filtered by lead time and/or admin level"""
        if not self.data_units:
            raise ValueError("Data units not found")
        if lead_time is not None and adm_level is not None:
            return list(
                filter(
                    lambda x: x.lead_time == lead_time and x.adm_level == adm_level,
                    self.data_units,
                )
            )
        elif lead_time is not None:
            return list(filter(lambda x: x.lead_time == lead_time, self.data_units))
        elif adm_level is not None:
            return list(filter(lambda x: x.adm_level == adm_level, self.data_units))
        else:
            return self.data_units

    def upsert_data_unit(self, data_unit: ClimateRegionDataUnit):
        """Add data unit; if it already exists, update it"""
        if not self.data_units:
            self.data_units = [data_unit]
        if hasattr(data_unit, "lead_time"):
            bdu = next(
                filter(
                    lambda x: x[1].climate_region_code == data_unit.climate_region_code
                    and x[1].lead_time == data_unit.lead_time,
                    enumerate(self.data_units),
                ),
                None,
            )
        else:
            bdu = next(
                filter(
                    lambda x: x[1].climate_region_code == data_unit.climate_region_code,
                    enumerate(self.data_units),
                ),
                None,
            )
        if not bdu:
            self.data_units.append(data_unit)
        else:
            self.data_units[bdu[0]] = data_unit

    def get_lead_times(self):
        """Return list of unique lead times"""
        return list(
            set([x.lead_time for x in self.data_units if hasattr(x, "lead_time")])
        )

    def get_climate_region_codes(self):
        """Return list of unique station codes"""
        return list(
            set(
                [
                    x.climate_region_code
                    for x in self.data_units
                    if hasattr(x, "climate_region_code")
                ]
            )
        )


class PipelineDataSets:
    """Collection of datasets used by the pipeline"""

    def __init__(
        self, country: str, settings: Settings, datetime: datetime = datetime.today()
    ):
        self.country = country

        self.rainfall_climateregion = ClimateRegionDataSet(
            country=self.country, timestamp=datetime
        )

        self.forecast_climateregion = ClimateRegionDataSet(
            country=self.country, timestamp=datetime
        )

        self.forecast_admin = AdminDataSet(
            country=self.country,
            timestamp=datetime,
            adm_levels=settings.get_country_setting(country, "admin-levels"),
        )

        self.threshold_climateregion = ClimateRegionDataSet(
            country=self.country, timestamp=datetime
        )
