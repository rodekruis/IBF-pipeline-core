
import pandas as pd
from datetime import datetime 
from azure.cosmos import CosmosClient
import geopandas as gpd
import os


COSMOS_URL=""
COSMOS_KEY=""
 
# Initialize the CosmosClient
client = CosmosClient(COSMOS_URL, credential=COSMOS_KEY, user_agent="ibf-flood-pipeline", user_agent_overwrite=True)

# Access a specific database
database_name = "drought-pipeline"
database = client.get_database_client(database_name)

country_admin_level = {"KEN": 1}#, "UGA": 2, "ZMB": 2}

for key,value in country_admin_level.items():
    country=key
    adm_level=value
    print(country,adm_level)

    df = pd.read_csv(f'data_updates/{country}_climate_region_district_mapping.csv')
    df['placeCode'] = df[f'ADM{adm_level}_PCODE'].astype(str)

    # Group the DataFrame by 'Climate_Region' and create a list of 'ParentPcode' for each group
    climate_region_parent_codes = df.groupby(['Climate_Region','Climate_Region_code'])['placeCode'].apply(list).to_dict()
    
    
    gd=gpd.read_file(f"data_updates/{country}.geojson")
    print(climate_region_parent_codes)
    print(gd.head())
    for climate_region, parent_codes in climate_region_parent_codes.items():
        gd_filtered=gd[gd[f'ADM{adm_level}_PCODE'].isin(parent_codes)]
        adm_dict = gd_filtered[['ADM1_PCODE', 'ADM2_PCODE','ADM3_PCODE']].to_dict(orient='list')
        adm_dict = {int(1): list(set(adm_dict['ADM1_PCODE'])), int(2): list(set(adm_dict['ADM2_PCODE'])), int(3): list(set(adm_dict['ADM3_PCODE']))}

        new_item={
            "country":country,
            "adm_level": adm_level,
            "climate_region":climate_region[0],
            "climate_region_code":climate_region[1],
            "pcodes":adm_dict,
            'timestamp': datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"),
            'id': f'{country}_{climate_region[1]}_{str(datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"))}',
        }
        container_name = "climate-region"
        container = database.get_container_client(container_name)
        container.upsert_item(new_item)