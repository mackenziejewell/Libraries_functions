
#////////////////////
#  shift_bylons  ///
#//////////////////
#---------------------------------------------------------------------
# Function to shift ERA5 data by longitudes. 
#---------------------------------------------------------------------
#///////////////////
#  sort_bylats  ///
#/////////////////
#---------------------------------------------------------------------
# Function to reverse ERA5 data by latitudes. 
#---------------------------------------------------------------------
#/////////////////
#  grab_ERA5  ///
#///////////////
#---------------------------------------------------------------------
# Function to grab ERA5 data and load into xarray dataframe.
#---------------------------------------------------------------------
#////////////////////
#  Download_ERA5 ///
#//////////////////
#---------------------------------------------------------------------
# Given datetime object, use cdsapi.Client() to download ERA5 variables
# This works for downloading hourly data from ERA5 on single levels. 
# For given date time, downloads from most recent hour before datetime.
#---------------------------------------------------------------------

#////////////////////
#  shift_bylons  ///
#//////////////////
#---------------------------------------------------------------------
# Function to shift ERA5 data by longitudes. 
#---------------------------------------------------------------------
# DEPENDENCIES
import xarray as xr
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def shift_bylons(ds, longitude_name = 'longitude', suppress_prints = True):

    """Function to shift ERA5 data by longitudes. 
    Given xarray-opened ECMWF dataset with longitude range [-180,180] range type, 
    shift to [0,360] range type. Only applies to data crossing the dateline, 
    with some positive and some negative longitudes. Helpful for errors in 
    plotting/calculating around discontinuous longitudes across pacific dateline.
    Shifted longitude coordinates will be saved in dataset as 'longitude'.

INPUT: 
- ds: xarray dataset with Nx1 longitude coordinate
- longitude_name: name of longitude coordinate in dataset 
                  (default: 'longitude')
- suppress_prints: whether or not to suppress operation descriptions 
                  (bool, default: True)

OUTPUT:
- new_ds: xarray dataset shifted about [0,360] range type 'longitude' coordinate

DEPENDENCIES:
import xarray as xr
import numpy as np
import numpy.ma as ma

Latest recorded update:
06-03-2022
    """

    # check if longitude range satisfies conditions
    if ds[longitude_name].values.min() < 0 and ds[longitude_name].values.max() > 0:
        if suppress_prints == False:
            print(f'Shifting longitudes from [{ds[longitude_name].values.min()}, {ds[longitude_name].values.max()}] range')

        # pull out longitude values as 'lons'
        lons=ds[longitude_name].values
        # add 360 to negative lons, shifting them [0,360]
        for ii in range(len(lons)):
            if lons[ii]<0:
                lons[ii]+=360   

        # assign new longitude coordinates from lons
        # variable will be 'longitude' (either replacing existing variable
        # or assigning new one if longitude_variable != 'longitude')
        new_lon_ds = ds.assign_coords(longitude = lons)
        new_ds = new_lon_ds.sortby('longitude')
        if suppress_prints == False:
            print(f">> 'longitude' variable range: [{new_ds['longitude'].values.min()}, {new_ds['longitude'].values.max()}]")

    else:
        if suppress_prints == False:
            print(f'Longitude range already satisfies conditions')
        # returning dataset as is
        new_ds = ds
            
    # return longitude-shifted data set
    return new_ds


#////////////////////
#  reverse_lats  ///
#//////////////////
#---------------------------------------------------------------------
# Function to reverse ERA5 data by latitudes. 
#---------------------------------------------------------------------
# DEPENDENCIES
import xarray as xr
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def reverse_lats(ds, latitude_name = 'latitude', suppress_prints = True):

    """Function to reverse ERA5 data by latitudes. 
    ERA5 data often given in descending order. If so, reverse to ascending order.

INPUT: 
- ds: xarray dataset with Nx1 longitude coordinate
- latitude_name: name of latitude coordinate in dataset 
                  (default: 'latitude')
- suppress_prints: whether or not to suppress operation descriptions 
                  (bool, default: True)

OUTPUT:
- new_ds: xarray dataset shifted to ascending order along 'latitude' coordinate

DEPENDENCIES:
import xarray as xr
import numpy as np
import numpy.ma as ma

Latest recorded update:
06-03-2022
    """

    # check if latitude is in descending order
    if ds[latitude_name][0].values>ds[latitude_name][-1].values:
        if suppress_prints == False:
            print(f'Reversing latitude from [{ds[latitude_name][0].values}, {ds[latitude_name][-1].values}] order')
        new_ds = ds.sortby(latitude_name)
        if suppress_prints == False:
            print(f">> '{latitude_name}' variable range: [{new_ds[latitude_name][0].values}, {new_ds[latitude_name][-1].values}]")

    else:
        if suppress_prints == False:
            print(f'Latitude range already satisfies conditions')
        # returning dataset as is
        new_ds = ds
            
    # return latitude-reversed data set
    return new_ds




#///////////////////
#  sort_bylats  ///
#/////////////////
#---------------------------------------------------------------------
# Function to reverse ERA5 data by latitudes. 
#---------------------------------------------------------------------
# DEPENDENCIES
import xarray as xr
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def sort_bylats(ds, latitude_name = 'latitude', suppress_prints = True):

    """Function to reverse ERA5 data by latitudes. 
    ERA5 data often given in descending order. If so, reverse to ascending order.

INPUT: 
- ds: xarray dataset with Nx1 longitude coordinate
- latitude_name: name of latitude coordinate in dataset 
                  (default: 'latitude')
- suppress_prints: whether or not to suppress operation descriptions 
                  (bool, default: True)

OUTPUT:
- new_ds: xarray dataset shifted to ascending order along 'latitude' coordinate

DEPENDENCIES:
import xarray as xr
import numpy as np
import numpy.ma as ma

Latest recorded update:
06-03-2022
    """

    # check if latitude is in descending order
    if ds[latitude_name][0].values>ds[latitude_name][-1].values:
        if suppress_prints == False:
            print(f'Reversing latitude from [{ds[latitude_name][0].values}, {ds[latitude_name][-1].values}] order')
        new_ds = ds.sortby(latitude_name)
        if suppress_prints == False:
            print(f">> '{latitude_name}' variable range: [{new_ds[latitude_name][0].values}, {new_ds[latitude_name][-1].values}]")

    else:
        if suppress_prints == False:
            print(f'Latitude range already satisfies conditions')
        # returning dataset as is
        new_ds = ds
            
    # return latitude-reversed data set
    return new_ds



#/////////////////
#  grab_ERA5  ///
#///////////////
#---------------------------------------------------------------------
# Function to grab ERA5 data and load into xarray dataframe.
#---------------------------------------------------------------------
# DEPENDENCIES
import xarray as xr
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def grab_ERA5(dt_obj, ERA5_path = '/Users/mackenziejewell/Data/ERA5.nc',
             lat_range = [65,90], lon_range = [100,300], crop_time = True):
    
    """Function to grab ERA5 data and load into xarray dataframe.

INPUT: 
- dt_obj: date of ERA5 data
- ERA5_path: path to ERA5 data (default: '/Users/mackenziejewell/Data/ERA5.nc')
- lat_range: [latmin, latmax] range to crop or None if should not crop (default = [65,90])
- lon_range: [lonmin, lonmax] range to crop or None if should not crop (default = [100,300])
- crop_time: whether or not to grab time of dt_obj from ds (True, default) or load in all time (False)

OUTPUT:
- new_ds: xarray dataset shifted to ascending order along 'latitude' coordinate

DEPENDENCIES:
import xarray as xr
import numpy as np
import numpy.ma as ma

Latest recorded update:
06-14-2022
    """
    
    # import data with xarray
    ds_unshifted = xr.open_dataset(ERA5_path)
    ds_unshifted.close() 
    
    # shift data to [0,360] longitude range (and ascending latitude range) if not already
    ERA5_ds = shift_bylons(ds_unshifted, longitude_name = 'longitude', suppress_prints=True)
    ERA5_ds = sort_bylats(ERA5_ds, latitude_name='latitude', suppress_prints=True)
    
    # grab single time
    if crop_time == True:
        ERA5_ds = ERA5_ds.sel(time = dt_obj)
    
    if lat_range != None:
        ERA5_ds = ERA5_ds.sel(latitude  = slice(lat_range[0],lat_range[-1]))
    if lon_range != None:
        ERA5_ds = ERA5_ds.sel(longitude = slice(lon_range[0],lon_range[-1])) 

    return ERA5_ds


#////////////////////
#  Download_ERA5 ///
#//////////////////
#---------------------------------------------------------------------
# Given datetime object, use cdsapi.Client() to download ERA5 variables
# This works for downloading hourly data from ERA5 on single levels. 
# For given date time, downloads from most recent hour before datetime.
#---------------------------------------------------------------------
# DEPENDENCIES
import datetime as dt
from datetime import datetime
import cdsapi
import xarray as xr
import numpy as np
import os
from urllib.request import urlopen
#---------------------------------------------------------------------

def Download_ERA5(dt_obj = None,
                  years = [2000],
                  months = np.arange(1,12+1),
                  days = np.arange(1,31+1),
                  hours = np.arange(1,23+1),
                  ERAprod = 'reanalysis-era5-single-levels',
                  ERAvars = ['10m_u_component_of_wind', '10m_v_component_of_wind','mean_sea_level_pressure'],
                  download_data = True, read_into_memory= False,
                  alt_dir = '',  SaveName = 'temp_ERAdownload.nc', 
                  extent = [90, -180, 60,180], quiet = False, allow_overwrites=False):
    
    """Function to download ERA5 data using cdsapi.Client(). Given datetime object or multiple dates information, automatically download ERA5 variables to file and/or memory. This works for downloading hourly data from ERA5 on single levels and may work for other data sets. When given single date time object, downloads from most recent hour before datetime.
       
       Referencing: https://towardsdatascience.com/read-era5-directly-into-memory-with-python-511a2740bba0

INPUT: 

IF DOWNLOADING SINGLE DATE (default, if dt_obj != None):
- dt_obj: datetime object of desired date to be found in ECMWF (default: None)
    (will download data from nearest hour before given time)
IF DOWNLOADING RANGE OF DATES (if dt_obj = None):
- years: list/array of years to download
- months: list/array of month numbers to download (default: np.arange(1,12+1), all months)
- days: list/array of day numbers to download (will automatically skip days not existing in month, for example day 30 of feb)
    (default: np.arange(1,31+1), all days)
- hours: list/array of hours to download (default: np.arange(1,23+1), all hours)
- ERAprod: ERA5 product to download (default: 'reanalysis-era5-single-levels')
- ERAvars: list of ERA5 variables to download 
            (default: ['10m_u_component_of_wind', '10m_v_component_of_wind','mean_sea_level_pressure'])
- download_data: bool, whether or not to download data (default: True)
- read_into_memory: bool, whether or not to read data into memory (default: False)
- alt_dir: directory to store downloaded .nc file (default: directory where code is run)
- SaveName: desired filename (default: 'temp_ERAdownload.nc')
- allow_overwrites: bool, whether or not to allow overwrites if file SaveName already exists
    (default: False)
- extent: data extent to download [N, W S, E] (default: [85, 170, 65, -110])
- quiet: option to quiet API client download information/automatic input
        (default: False)

OUTPUT:
if download_data == True and read_into_memory == False:
    - FullFileName: Name and directory where file was downloaded
elif download_data == False and read_into_memory == True:
    - ds: xarray dataset of downloaded data
elif download_data == True and read_into_memory == True:
    - FullFileName: Name and directory where file was downloaded
    - ds: xarray dataset of downloaded data

DEPENDENCIES:
import datetime as dt
from datetime import datetime
import cdsapi
import xarray as xr
import numpy as np
import os
from urllib.request import urlopen

Latest recorded update:
06-23-2023
    """

    # check whether file already exists if downloading
    allow_download = True
    if download_data == True:
        FullFileName = alt_dir+SaveName
        if os.path.isfile(FullFileName):
            if allow_overwrites:
                if quiet == False:
                    print(f' >>> {SaveName} already exists in {alt_dir}. Will overwrite file.')
            else:
                print(f' >>> {SaveName} already exists in {alt_dir}. Will not overwrite file.')
                allow_download = False
    
    # check for dt_obj. If not None, use this download date. 
    if dt_obj != None:
        
        # grab date data from datetime object
        #------------------------------------
        Year_string = dt_obj.strftime('%Y')
        Month_string = dt_obj.strftime('%m')
        Day_string = dt_obj.strftime('%d')
        Time_string = dt_obj.strftime('%H')+':00'
        NearestDate = datetime.strptime(Year_string+Month_string+Day_string+dt_obj.strftime('%H'), '%Y%m%d%H')
        
        if quiet == False:
            print(f' >>> download: {NearestDate}')

    # otherwise download across multiple dates as below
    else:  
        Year_string = [str(year) for year in years]
        Month_string = [str(month).zfill(2) for month in months]
        Day_string = [str(day).zfill(2) for day in days]
        Time_string = [f'{str(hours).zfill(2)}:00' for hours in hours]
        
        if quiet == False:
            print(f' >>> download across dates:')
            print(f'     - years: {Year_string}')
            print(f'     - months: {Month_string}')
            print(f'     - days: {Day_string}')
            print(f'     - hours: {Time_string}')


    # make download with cdsapi.Client()
    # if exception occurs regarding Missing/incomplete configuration file,
    # print instructions to set one of these up
    #-----------------------------------
    error = None
    try:
        c = cdsapi.Client(quiet = quiet)
    except Exception as e: 
        error = e
        print(e)
    if error != None:
        if 'Missing/incomplete configuration file' in str(error):
            print(f'   make sure CDS API configuration file has been set up on computer.')
            print(f'   Reference {"https://cds.climate.copernicus.eu/api-how-to"} for set-up information.')
            print(f'   {"https://bookdown.org/huckley/microclimate_users_guide/era5.html"} also has more detailed instructions.')
            
    # api parameters 
    params = {'product_type': 'reanalysis',
              'variable': ERAvars,
              'year': Year_string,
              'month': Month_string,
              'day': Day_string,
              'time': Time_string,
              'area': extent,
              'format': 'netcdf'}
    
    # retrieve the path to the file
    #------------------------------
    data = c.retrieve(ERAprod, params)

    # download the file, if specified
    #--------------------------------
    if download_data and allow_download:
        data.download(FullFileName)
        if quiet == False:
            print(f' >>> save data to: {FullFileName}')
    
    # read data into memory, if specified
    #------------------------------------
    if read_into_memory:
        with urlopen(data.location) as f:
            ds = xr.open_dataset(f.read())
    
    # return downloaded data path, or dataset in memory, or both.
    #------------------------------------------------------------
    if download_data == True and read_into_memory == False:
        return FullFileName
    elif download_data == False and read_into_memory == True:
        return ds
    elif download_data == True and read_into_memory == True:
        return FullFileName, ds
    else: 
        print('At least one of download_data and read_into_memory must be True.')