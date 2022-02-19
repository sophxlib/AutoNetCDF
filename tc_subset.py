'''
Hurricane Data Output
Purpose: Create a netCDF output from IBTracks NOAA MIRS satellite data
User input: Name of storm and year of storm 
Note to user: Change lines 48 to put MIRS directory path
Comment out lines 71,72 if you are testing or running production
'''
import pandas as pd 
import pytz
import os
from shapely import wkt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import xarray as xr


def LonTo360(dlon):
    # Convert longitudes to 0-360 deg
    dlon = ((360 + (dlon % 360)) % 360)
    return dlon
    
def get_by_name(name: str, year: int, ibtrack_file: str):
    data = pd.read_csv(ibtrack_file, low_memory=False)  
    data = data.iloc[1: , :] # remove the row of units
    year_date_of = pd.to_datetime(f'{year}')
    year_date_after = pd.to_datetime(f'{year+1}')
    data = data[data['NAME'] == name]
    data['ISO_TIME'] = pd.to_datetime(data['ISO_TIME'])
    mask = (data['ISO_TIME'] >= year_date_of) & (data['ISO_TIME'] < year_date_after)
    data = data[mask]
    extract_vars = ['NAME', 'ISO_TIME', 'WMO_WIND', 'WMO_PRES', 'LAT', 'LON', 
                    'USA_R34_NE', 'USA_R34_NW', 'USA_R34_SE', 'USA_R34_SW', 'USA_R50_NE', 'USA_R50_NW', 'USA_R50_SE', 'USA_R50_SW', 'USA_R64_NE', 'USA_R64_NW', 'USA_R64_SE', 'USA_R64_SW', 
                   'REUNION_R34_NE', 'REUNION_R34_NW', 'REUNION_R34_SE', 'REUNION_R34_SW', 'REUNION_R50_NE', 'REUNION_R50_NW', 'REUNION_R50_SE', 'REUNION_R50_SW', 'REUNION_R64_NE', 'REUNION_R64_NW', 'REUNION_R64_SE', 'REUNION_R64_SW',
                   'BOM_R34_NE', 'BOM_R34_SE', 'BOM_R34_NW', 'BOM_R34_SW', 'BOM_R50_NE', 'BOM_R50_SE', 'BOM_R50_NW', 'BOM_R50_SW', 'BOM_R64_NE', 'BOM_R64_SE', 'BOM_R64_NW', 'BOM_R64_SW']
    data = data[extract_vars]
    if year != 2021: 
        data = data[data['WMO_WIND'] != ' ']
        data = data[data['WMO_PRES'] != ' ']
    data['LON_180'] = data['LON']
    data['LON']=data['LON'].astype(float).apply(LonTo360)
    return data

### user settings
storm_name = "IDA"
storm_year = 2021
output_filename = f"{storm_name}_{storm_year}_all_data.nc"

dir_path = '/Users/sophiahu/Documents/MIRS_DATA/' # change to location of MIRS data
ibt_file = 'IBTrACS.ALL.v04r00.nc'
ibt_file_csv = 'ibtracs.ALL.list.v04r00.csv'

snd_vars_to_keep = ["Player", "Plevel", "PTemp", "PVapor", "PClw", "PRain", "PGraupel"] # unique from IMG

img_vars_to_remove = ["Atm_type", "ChanSel", "SWP", "IWP", "Snow",
                      "SWE", "SnowGS", "SIce", "SIce_MY", "SIce_FY", "SFR",
                     "CldTop", "CldBase", "CldThick", "PrecipType", "RFlag", "SurfM",
                     "WindSp", "WindDir", "WindU", "WindV", "Prob_SF", "quality_information"]

################
## main code
################

mirs_files = os.listdir(dir_path)
result = get_by_name(storm_name, storm_year, ibt_file_csv)


npts = result.shape[0] # read result DataFrame for storm name and year
mirs_contains_tc = [] # list of MIRS files that contain TC 


#for ipt in range(npts):    # turn this on in production; (loop through all points of storm's lifecyle)
for ipt in [42,43]:         # turn this off in production; testing only (loop through a few points of storm's lifecycle - here 42, 43 are for 08/29/21 for IDA
    
    t_val = result['ISO_TIME'].iloc[ipt] 
    lat_val = result['LAT'].iloc[ipt] 
    lon_val = result['LON_180'].iloc[ipt] # LON is from 0-360 but the satellite data uses -180 to 180!

    print(f"Working on TC time: {t_val}")

    t_val_12ahead = t_val + pd.to_timedelta(12, unit='h')
    t_val_12behind = t_val - pd.to_timedelta(12, unit='h')

    t_val_12ahead_utc = pytz.utc.localize(t_val_12ahead)
    t_val_12behind_utc = pytz.utc.localize(t_val_12behind)

    print(f"Checking files within date range: {t_val_12behind} and {t_val_12ahead}")
        
    for file in mirs_files: 
        ds = xr.open_dataset(dir_path + file)
        time_coverage_start = pd.to_datetime(ds.attrs['time_coverage_start'])
        time_coverage_end = pd.to_datetime(ds.attrs['time_coverage_end'])

        if (t_val_12behind_utc <= time_coverage_start) and (t_val_12ahead_utc >= time_coverage_end):
            polygon= ds.attrs['geospatial_bounds']
            polygon_val = wkt.loads(polygon)
            point = Point(float(lon_val), float(lat_val))
            
            if point.within(polygon_val): 
                if ds.attrs["geospatial_first_scanline_first_fov_lon"] < 0:  # keep storms East of the dateline for IDA; may need to alter for other storms (to remove satellite granulars that circle the globe)
                    mirs_contains_tc.append(file)
                    
img_files_storm_final = [i for i in mirs_contains_tc if i.startswith('NPR-MIRS-IMG')]
snd_files_storm_final = [i for i in mirs_contains_tc if i.startswith('NPR-MIRS-SND')]                    
                    
# Concatenate the MIRS files
img_files_storm_final_one_list = [item for sublist in img_files_storm_final for item in sublist]
ds_img_list = []

for file in img_files_storm_final_one_list:
    ds_mirs = xr.open_dataset(dir_path + file)
    ds_img_list.append(ds_mirs)

ds_img_merged = xr.concat(ds_img_list, dim='Scanline')

snd_files_storm_final_one_list = [item for sublist in snd_files_storm_final for item in sublist]
ds_snd_list = []

for file in snd_files_storm_final_one_list:
    ds_mirs = xr.open_dataset(dir_path + file)
    ds_snd_list.append(ds_mirs)

ds_snd_merged = xr.concat(ds_snd_list, dim='Scanline')


# Keep certain SND variables and remove certain IMG variables
ds_snd_merged_keep_vars = ds_snd_merged[snd_vars_to_keep]
ds_img_merged_keep_vars = ds_img_merged.drop(img_vars_to_remove)

# Read the IBtracks data for the storm
ds_ibt = xr.open_dataset(ibt_file)
storm_name_bytes = bytes(storm_name, 'UTF-8')
ds_storm_all = ds_ibt.where(ds_ibt.name==storm_name_bytes, drop=True)
ds_storm = ds_storm_all.where(ds_storm_all.season==float(storm_year), drop=True)

# Merge the IBtracks and MIRS data into one file
ds_merged_all = xr.merge([ds_img_merged_keep_vars, ds_snd_merged_keep_vars, ds_storm])


# Add attributes
ds_merged_all.attrs["TC_name"] = storm_name
ds_merged_all.attrs["TC_time_start"] = bytes.decode( ds_storm["iso_time"][:,0].item() )
# ds_merged_all.attrs["TC_time_end"] =  bytes.decode( ds_storm["iso_time"][:,-1].item() ) # len_iso_time = len(i)
ds_merged_all.attrs["TC_minimum_lat"] = round(float(ds_storm["lat"].min()),2)
ds_merged_all.attrs["TC_minimum_lon"] = round(float(ds_storm["lon"].min()),2)
ds_merged_all.attrs["TC_maximum_lat"] = round(float(ds_storm["lat"].max()),2)
ds_merged_all.attrs["TC_maximum_lon"] = round(float(ds_storm["lon"].max()),2)

ds_merged_all.to_netcdf(output_filename)