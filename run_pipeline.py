import logging

from ibfpipeline.riverflood.pipeline import Pipeline as RiverFloodPipeline
from ibfpipeline.drought.pipeline import Pipeline as DroughtPipeline
from ibfpipeline.core.secrets import Secrets
from ibfpipeline.core.settings import Settings
import click


@click.command()
@click.option("--hazard", help="hazard name", default="riverflood")
@click.option("--country", help="country ISO3", default="UGA")
@click.option("--prepare", help="prepare discharge data", default=False, is_flag=True)
@click.option("--forecast", help="forecast floods", default=False, is_flag=True)
@click.option("--send", help="send to IBF", default=False, is_flag=True)
@click.option(
    "--debug",
    help="debug mode: process only one ensemble member from yesterday",
    default=False,
    is_flag=True,
)
def run_pipeline(hazard, country, prepare, forecast, send, debug):
    country = country.upper()
    try:
        if hazard.lower() == "riverflood":
            pipe = RiverFloodPipeline(
                country=country,
                settings=Settings("config/riverflood/config.yaml"),
                secrets=Secrets(".env"),
            )
        elif hazard.lower() == "drought":
            pipe = DroughtPipeline(
                country=country,
                settings=Settings("config/drought/config.yaml"),
                secrets=Secrets(".env"),
            )
        else:
            raise ValueError(f"Hazard {hazard} not supported.")
    except FileNotFoundError as e:
        logging.warning(f"Necessary dataset missing: {e}, skipping country {country}")
        return

    pipe.run_pipeline(
        prepare=prepare,
        forecast=forecast,
        send=send,
        debug=debug,
    )


if __name__ == "__main__":
    run_pipeline()
