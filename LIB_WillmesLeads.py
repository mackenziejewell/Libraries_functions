import os
import shutil
import numpy as np, numpy.ma as ma
import xarray as xr

def grab_WillmesLeads(dt_obj, main_path='/Users/mackenziejewell/Data/willmesetal_leads_daily/',
                      doi = '10.1594/PANGAEA.854411', 
                      lon_range=[0,360], lat_range=[60,90], return_vars = ['lon', 'lat', 'ice_data', 'ds']):
    
    """Function to unzip, open, and return Willmes lead data stored locally on computer. 
    From: https://doi.pangaea.de/10.1594/PANGAEA.854411 or https://doi.pangaea.de/10.1594/PANGAEA.955561
    Files should be stored by year in folders labeled as ArcLeads_{} where {} is the year
    (in YYYY format). Each annual data folder should be in given main_path. Associated geo grid
    for lead data should be saved as 'latlonmap.nc' in main_path.
    

INPUT: 
- dt_obj: given date (as datetime object)
- main_path: path where individual ArcLeads folders are stored and latlonmap.nc file is stored (for version 854411).
             path where individual ArcLeads files are stored (for version 955561).
- doi: doi associated with version of ARCleads product (files stored differently).  
    - doi = '10.1594/PANGAEA.854411' for initial 2015 version where daily data are stored 
            in annual folders and geo coordinates stored separately (latlonmap.nc)
            (https://doi.org/10.1594/PANGAEA.854411)
    - doi = '10.1594/PANGAEA.955561' for 2023 version where data and geocoords are stored in annual files.
            (https://doi.org/10.1594/PANGAEA.955561)
- lon_range: range to crop data as [lonmin, lonmax] (default: [0,360] i.e. no cropping)
- lat_range: range to crop data as [latmin, latmax] (default: [60,90])
- return_vars: list of variables/attributes to return in specified order
    Can include any or all OUTPUT variables in any order: 
    default: ['lon', 'lat', 'freq', 'ds']
    - 'lon' : (M x N array of longitudes)
    - 'lat' : (M x N array of latitudes)
    - 'ice_data' : (M x N array of ice data classifications) (if not data available on this date, returns empty list)
    - 'ds' : xarray dataset (if not data available on this date, returns empty list)


OUTPUT:
- return_data: list of all data variables corresponding to names specified in return_vars
    for default return_vars specified above, read all output into list of variables for example as [lat, lon, ice_data, ds]

example:
[lat, lon, ice_data, ds] = grab_WillmesLeads(return_vars = ['lon', 'lat', 'ice_data', 'ds'])

DEPENDENCIES:
from os.path import exists
import shutil
import numpy as np, numpy.ma as ma
import xarray as xr

Latest recorded update:
06-20-2023
   """

    
    lat = []
    lon = []
    
    assert type(doi) == str, f"doi should be string, not {type(doi)}"
    assert doi in ['10.1594/PANGAEA.854411', '10.1594/PANGAEA.955561'], f"unrecognized doi, should be '10.1594/PANGAEA.854411' or '10.1594/PANGAEA.955561' "
    
    # store all output in dictionary
    vars_dict = {}
    
    if str(doi) == '10.1594/PANGAEA.854411':
        assert os.path.exists(main_path+'latlonmap.nc'), f"latlonmap.nc not found in {main_path}"

        # open geographic data coordinates from willmes lead data
        # 'latlonmap.nc' file
        #--------------------------------------------------------
        ds_map = xr.open_dataset(main_path+'latlonmap.nc')
        lat = ds_map.latitude.values
        lon = ds_map.longitude.values
        lon[lon<=0]+=360 # shift to (0,360) range
   
        # open folder associated with year of given date
        #-----------------------------------------------
        year = dt_obj.year
        year_folder = f'ArcLeads_{str(year)}/'
        file = 'ArcLeads_{}{}{}.nc'.format(dt_obj.strftime('%Y'),dt_obj.strftime('%m'),dt_obj.strftime('%d'))
        file_with_path = main_path+year_folder+file

    
        # check if file needs to be unzipped
        # if so, use shutil to unzip
        #-------------------------------------------------

        # check if file already exists as unzipped data.
        unzipped_file_exists = os.path.exists(file_with_path)

        # Full path of the archive file
        zipped_filename = file_with_path+'.zip'
        zipped_file_exists = os.path.exists(zipped_filename)


        # if file does not exist as zipped or unzipped file, break and return empty arrays
        if unzipped_file_exists == False and zipped_file_exists == False:
            print(f'>>> {file_with_path} does not exist')
            # read data into dict
            vars_dict['ds'] = []
            vars_dict['ice_data'] = []

        # if file exists, zipped or unzipped, keep running to extract data
        else:

            # if file has not yet been unzipped, unzip it
            if unzipped_file_exists == False and zipped_file_exists == True:
                print(f'>> unzip {file}')
                # Target directory
                extract_dir = main_path+year_folder
                # Unpack the archive file
                shutil.unpack_archive(zipped_filename, extract_dir, "zip")

            # open file, which should be unzipped if it wasn't already
            #---------------------------------------------------------
            ds = xr.open_dataset(file_with_path)
            leads = ds.leadMap.values  

    
    elif str(doi) == '10.1594/PANGAEA.955561':
        
        # construct filename for file containing dt_obj
        # months Jan - Apr of year YYYY are stored in file YYYY-1/YYYY
        # but listed with times YYYY-1 in the netcdf file
        dsyear = int(dt_obj.strftime('%Y'))
        if int(dt_obj.strftime('%m')) < 6:
            dsyear-=1
        file = f"{dsyear}{str(dsyear+1)[2:]}_ArcLeads.nc"
        file_with_path = main_path+file
        
        # check if file exists
        file_exists = os.path.exists(file_with_path)
        
        # if file does not exist, save empty lists
        if file_exists == False:
            print(f'>>> {file_with_path} does not exist')
            # read data into dict
            vars_dict['ds'] = []
            vars_dict['ice_data'] = []
           
        # if file exists, open file
        else:
            ds = xr.open_dataset(file_with_path)
            dstime = int( str(dsyear) + dt_obj.strftime('%m%d') )
            leads = ds.leadmap[np.where(ds.time==dstime)][0,:,:].values  
            lat = ds.lat.values
            lon = ds.lon.values
            lon[lon<=0]+=360 # shift to (0,360) range
            
    # check that lat data exists, if so, crop:
    if len(lat) > 0:
        
        # crop map as given by lon_range and lat_range
        #---------------------------------------------
        rows_to_delete = []
        columns_to_delete = []

        # check for values within ranges row by row
        for ii in range(np.shape(lon)[0]):
            # check how many indices along iith row are within lon_range and lat_range
            good_lons_ii = np.sum((lon[ii,:]>=lon_range[0]) & (lon[ii,:]<=lon_range[1]))
            good_lats_ii = np.sum((lat[ii,:]>=lat_range[0]) & (lat[ii,:]<=lat_range[1]))
            # if no lons or no lats are within range, delete row
            if good_lons_ii == 0 or good_lats_ii == 0:
                rows_to_delete.append(ii)

        # check for values within ranges column by column
        for jj in range(np.shape(lon)[1]):
            # check how many indices along iith row are within lon_range and lat_range
            good_lons_jj = np.sum((lon[:,jj]>=lon_range[0]) & (lon[:,jj]<=lon_range[1]))
            good_lats_jj = np.sum((lat[:,jj]>=lat_range[0]) & (lat[:,jj]<=lat_range[1]))
            # if no lons or no lats are within range, delete column
            if good_lons_jj == 0 or good_lats_jj == 0:
                columns_to_delete.append(jj)

        lats = np.delete(lat, columns_to_delete, 1)
        lats = np.delete(lats, rows_to_delete, 0)
        lons = np.delete(lon, columns_to_delete, 1)
        lons = np.delete(lons, rows_to_delete, 0)

        # read geo data into dict
        vars_dict['lat'] = lats
        vars_dict['lon'] = lons
        
        # delete values outside lon_range and lat_range to match geo data
        ice_data = np.delete(leads, columns_to_delete, 1)
        ice_data = np.delete(ice_data, rows_to_delete, 0)

        # read data into dict
        vars_dict['ds'] = ds
        vars_dict['ice_data'] = ice_data
    
    else:
        vars_dict['lat'] = []
        vars_dict['lon'] = []
    
    
    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    if len(return_data) == 1:
        return_data = return_data[0]
    
    return return_data
           