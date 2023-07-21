#//////////////////////////
#  convert_PPD_vectors ///
#////////////////////////
#---------------------------------------------------------------------
# Convert u/v (polar EASE grid) velocity components from NSIDC PPD
# data to eastward and northward vector components.
#---------------------------------------------------------------------
#////////////////////
#  crop_PPD_data ///
#//////////////////
#---------------------------------------------------------------------
# Crop NSIDC Polar Pathfinder lats, lons, u, v, to within given lat/lon range
#---------------------------------------------------------------------
#/////////////////////
#  grab_ice_Drift ///
#///////////////////
#---------------------------------------------------------------------
# Import NSIDC Polar Pathfinder lats, lons, u, v cropped to within given lat/lon range.
#---------------------------------------------------------------------
#//////////////////////////
#  grab_icedrift_range ///
#////////////////////////
#---------------------------------------------------------------------
# Import NSIDC Polar Pathfinder data over desired time range within year.
#---------------------------------------------------------------------




#//////////////////////////
#  convert_PPD_vectors ///
#////////////////////////
#---------------------------------------------------------------------
# Convert u/v (polar EASE grid) velocity components from NSIDC PPD
# data to eastward and northward vector components.
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
#---------------------------------------------------------------------
def convert_PPD_vectors(lon = [], u_EASE = [], v_EASE = [], fill_to_nan = False):
    
    """Convert u/v (polar EASE grid) velocity components from NSIDC Polar Pathfinder data to 
    eastward and northward vector components. ref: https://nsidc.org/support/how/how-convert-horizontal-and-vertical-components-east-and-north.
    
INPUT:
- u_EASE: M x N grid of along-x velocity component on EASE grid (cm/s)
- v_EASE: M x N grid of along-y velocity component on EASE grid (cm/s)
    (u_EASE: toward the right on the grid)
    (v_EASE: upward (toward the top) on the grid)
- lon: M x N grid of longitude grid (0 to 360) associated with u, v
- fill_to_nan: specify whether or not to replace fill value in masked arrays 
    with nans in calculation so fill number is not treated as actual data value. 
    (default: False if no changes applied to data)
    (else set to True and routine will determine masked array fill value and replace with nans).

OUTPUT:
- u:  M x N grid of east component of velocity (cm/s)
- v:  M x N grid of north component of velocity (cm/s)

DEPENDENCIES:
import numpy as np

Latest recorded update:
06-28-2023
    """
                
    if fill_to_nan:
        u_fill = u_EASE.fill_value
        u_EASE = u_EASE.data
        u_EASE[u_EASE == u_fill] = np.nan
        
        v_fill = v_EASE.fill_value
        v_EASE = v_EASE.data
        v_EASE[v_EASE == fill_to_nan] = np.nan
        
    # convert EASE grid vector components to northward, eastward
    #------------------------------------------------------------------
    u = u_EASE * np.cos(lon/180*np.pi)  +  v_EASE * np.sin(lon/180*np.pi)
    v = -u_EASE * np.sin(lon/180*np.pi)  +  v_EASE * np.cos(lon/180*np.pi)
                
    return u, v



#////////////////////
#  crop_PPD_data ///
#//////////////////
#---------------------------------------------------------------------
# Crop NSIDC Polar Pathfinder lat, lon, u, v, to within given lat/lon range
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def crop_PPD_data(lon = [], lat = [], VARS = [], lon_range=[195, 235], lat_range=[65,78]):
    
    """Crop NSIDC Polar Pathfinder lat, lon, u, v, to within given lat/lon range.
    
INPUT:
- lon: M x N longitude grid (0 to 360) associated with u, v
- lat: M x N latitude grid associated with u, v
- VARS: list of variables to crop (default: empty list), each as a M x N grid
- lon_range:longitude range for cropping (defined 0 to 360)
    (default: [195, 235])
- lat_range:latitude range for cropping
    (default: [65,78])

OUTPUT:
- lon, lat, VARS  --> cropped input grids
    
DEPENDENCIES:
import numpy as np
import numpy.ma as ma

Latest recorded update:
07-21-2023
    """
    
    # grab min and max lat/lon values
    #--------------------------------------------------
    lonmin, lonmax = lon_range[0], lon_range[1]
    latmin, latmax = lat_range[0], lat_range[1]
    
    # run through lat/lon arrays and crop repeatedly 
    # until array sizes no longer change
    #-----------------------------------------------
    array_size_change = 10
    while array_size_change>0:
        # determine current size of arrays
        num_rows_i = np.shape(lat)[0]
        num_cols_i = np.shape(lat)[1]
        # run through each row in lat/lon and check 
        # if any coords are within desired range, keep row
        #--------------------------------------------------
        rows_to_delete = []
        for ii in range(np.shape(lat)[0]):
            # check whether all lat/lon in row are within ranges
            row_check = 0
            for jj in range(np.shape(lat)[1]):
                # add 1 to row_check and break current loop if any coords are within ranges
                if lat[ii][jj] < latmax and lat[ii][jj] > latmin and lon[ii][jj] < lonmax and lon[ii][jj] > lonmin:
                    row_check+=1
                    break
            # if no coords values within desired range were found in row,
            # add row index to rows_to_delete
            if row_check == 0:
                rows_to_delete.append(ii)

        # crop lat, lon, u, v to new range
        #------------------------------------
        lat = np.delete(lat, rows_to_delete, axis=0)
        lon = np.delete(lon, rows_to_delete, axis=0)
        
        for vv in range(len(VARS)):
            # VARS[vv] can have dims (time, lat, lon) or just (lat, lon)
            # either way, find axis of lat to crop
            AX = len(VARS[vv].shape)-2
            VARS[vv] = np.delete(VARS[vv], rows_to_delete, axis=AX)

        # run through each column in lat/lon and check 
        # if any coords are within desired range, keep column
        #----------------------------------------------------
        columns_to_delete = []
        for jj in range(np.shape(lat)[1]):
            # check whether all lat/lon in column are within ranges
            column_check = 0
            for ii in range(np.shape(lat)[0]):
                # add 1 to column_check and break current loop if any
                # coords are within ranges
                if lat[ii][jj] < latmax and lat[ii][jj] > latmin and lon[ii][jj] < lonmax and lon[ii][jj] > lonmin:
                    column_check+=1
                    break
            # if no coords values within desired range were found in column, add column index to columns_to_delete
            if column_check == 0:
                columns_to_delete.append(jj)

        # crop lat, lon, u, v to new range
        #------------------------------------
        lat = np.delete(lat, columns_to_delete, axis=1)
        lon = np.delete(lon, columns_to_delete, axis=1)
        for vv in range(len(VARS)):
            # VARS[vv] can have dims (time, lat, lon) or just (lat, lon)
            # either way, find axis of lon to crop
            AX = len(VARS[vv].shape)-1
            VARS[vv] = np.delete(VARS[vv], columns_to_delete, axis=AX)

        # determine change in array sizes to see if it is still cropping
        array_size_change = (num_rows_i-np.shape(lat)[0])+(num_cols_i-np.shape(lat)[1])

    # np.delete messes with mask, so remask
    for vv in range(len(VARS)):
        VARS[vv] = ma.masked_where(VARS[vv] == -9999, VARS[vv])
    
    return lon, lat, VARS



#/////////////////////
#  grab_ice_Drift ///
#///////////////////
#---------------------------------------------------------------------
# Import NSIDC Polar Pathfinder lats, lons, u, v cropped to within given lat/lon range.
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
import xarray as xr
from datetime import datetime
#---------------------------------------------------------------------
def grab_ice_Drift(date = datetime(year = 2000, month = 1, day = 1),
                   PPD_drift_path = '/Volumes/Jewell_EasyStore/NSIDC-0116_PPdrift/', 
                   PPD_filename = 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc', 
                   return_vars = ['lon', 'lat', 'u', 'v', 'xx', 'yy', 'u_EASE', 'v_EASE', 'proj', 'ds'], 
                   lat_range = [0, 90], lon_range = [0, 360]):

    """Import NSIDC Polar Pathfinder (sea ice drift NSIDC-0116, doi:10.5067/INAWUWO7QH7B) lats, lons, u, v cropped to within given lat/lon range.

INPUT:
- date: desired date (datetime object)
- PPD_drift_path: directory where PPD files are locally stored.
- PPD_filename: naming convention for PPD files (default: 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc' where {} will be replaced with year of dt_obj)
- return_vars: variables/attributes to return in specified order. (list)
    Can include any or all OUTPUT variables in any order: 
    default: ['lon', 'lat', 'u', 'v', 'xx', 'yy', 'u_EASE', 'v_EASE', 'proj', 'ds']
- lat_range:latitude range for cropping
    (default: [0, 90])
- lon_range:longitude range for cropping (defined 0 to 360)
    (default: [0, 360])

OUTPUT:
List of any or all of variables specified in return_vars:
- lon: M x N longitude grid (0 to 360) associated with u, v
- lat: M x N latitude grid associated with u, v
- u: M x N grid of eastward vector components of ice drift
- v: M x N grid of northward vector components of ice drift
- xx: M x N grid of x values of EASEgrid, corresponding to u_EASE, v_EASE
- yy: M x N grid of y values of EASEgrid, corresponding to u_EASE, v_EASE
- u_EASE: M x N grid of along-x component of ice drift
- v_EASE: M x N grid of along-y component of ice drift
- proj: cartopy projection from PP drift data projection info
- ds: xarray data frame containing data from year of date

DEPENDENCIES:
import numpy as np
import numpy.ma as ma
import xarray as xr
from datetime import datetime

Latest recorded update:
06-28-2022
    """


    # store all output in dictionary
    vars_dict = {}

    # incomplete year of data in 1978
    if date.year == 1978:
        PPD_filename = 'icemotion_daily_nh_25km_19781101_19781231_v4.1.nc'

    # open file corresponding to date's year
    PP_mainpath = PPD_drift_path+PPD_filename
    PP_file = PP_mainpath.format(date.year,date.year)
    ds = xr.open_dataset(PP_file)
    ds.close()
    vars_dict['ds'] = ds

    # cartopy projection for xx, yy, u_EASE, v_EASE
    vars_dict['proj'] = grab_projinfo_PPdrift(ds)

    # find date match in data
    date_index = False
    for tt, time in enumerate(ds.time.values):
        if time.year == date.year:
            if time.month == date.month:
                if time.day == date.day:
                    date_index = tt               

    # if date match found, select data from specific day                
    if type(date_index) == bool:
        raise Exception(f"No date match found for {date.date()} in {PP_file}")
    else:
        ds_date = ds.sel(time = ds.time[date_index])

    # projected drift components
    u_EASE = ds_date.u.values
    v_EASE = ds_date.v.values

    # projected coords
    xx, yy = np.meshgrid(ds_date.x.values, ds_date.y.values)

    # lat/lon coords
    lat = ds_date.latitude.values
    lon = ds_date.longitude.values
    lon[lon<0]+=360

    # crop data to lat/long range
    lon, lat, [xx, yy, u_EASE, v_EASE] = crop_PPD_data(lon = lon, lat = lat, VARS = [xx, yy, u_EASE, v_EASE], lon_range=lon_range, lat_range=lat_range)

    # add to dict
    vars_dict['lon'] = lon
    vars_dict['lat'] = lat
    vars_dict['xx'] = xx
    vars_dict['yy'] = yy
    vars_dict['u_EASE'] = u_EASE
    vars_dict['v_EASE'] = v_EASE

    # if converted u, v desired, crop them
    if 'u' in return_vars or 'v' in return_vars:
        u, v = convert_PPD_vectors(lon = lon, u_EASE = u_EASE, v_EASE = v_EASE, fill_to_nan = True)
        vars_dict['u'] = u
        vars_dict['v'] = v

    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    if len(return_data) == 1:
        return_data = return_data[0]
    return return_data
            

#//////////////////////////
#  grab_icedrift_range ///
#////////////////////////
#---------------------------------------------------------------------
# Import NSIDC Polar Pathfinder data over desired time range within year.
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
import xarray as xr
from datetime import datetime
import pandas as pd
#---------------------------------------------------------------------
def grab_icedrift_range(start_date = datetime(year = 2000, month = 1, day = 1),
                        end_date = datetime(year = 2000, month = 12, day = 31),
                        PPD_drift_path = '/Volumes/Jewell_EasyStore/NSIDC-0116_PPdrift/', 
                        PPD_filename = 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc', 
                        return_vars = ['lon', 'lat', 'u', 'v', 'xx', 'yy', 'u_EASE', 'v_EASE', 'proj', 'ds'], 
                        lat_range = [0, 90], lon_range = [0, 360]):

    """Import NSIDC Polar Pathfinder (sea ice drift NSIDC-0116, doi:10.5067/INAWUWO7QH7B) data over time range.
        Time range must be within same year. Return xarray ds, or lats, lons, u, v cropped to within given lat/lon range.

INPUT:
- start_date: datetime object, initial date of desired period.
- end_date: datetime object, final date of desired period. Must be within same year as start_date.
- PPD_drift_path: directory where PPD files are locally stored.
- PPD_filename: naming convention for PPD files (default: 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc' where {} will be replaced with year of dt_obj)
- return_vars: variables/attributes to return in specified order. (list)
    Can include any or all OUTPUT variables in any order: 
    default: ['lon', 'lat', 'u', 'v', 'xx', 'yy', 'u_EASE', 'v_EASE', 'proj', 'ds']
- lat_range:latitude range for cropping
    (default: [0, 90])
- lon_range:longitude range for cropping (defined 0 to 360)
    (default: [0, 360])

OUTPUT:
List of any or all of variables specified in return_vars:
- lon: M x N longitude grid (0 to 360) associated with u, v
- lat: M x N latitude grid associated with u, v
- u: M x N grid of eastward vector components of ice drift
- v: M x N grid of northward vector components of ice drift
- xx: M x N grid of x values of EASEgrid, corresponding to u_EASE, v_EASE
- yy: M x N grid of y values of EASEgrid, corresponding to u_EASE, v_EASE
- u_EASE: M x N grid of along-x component of ice drift
- v_EASE: M x N grid of along-y component of ice drift
- proj: cartopy projection from PP drift data projection info
- ds: xarray data frame containing data from year of date

DEPENDENCIES:
import numpy as np
import numpy.ma as ma
import xarray as xr
from datetime import datetime
import pandas as pd

Latest recorded update:
07-21-2023
    """


    assert start_date.year == end_date.year, f"start year ({start_date.year }) and end year ({end_date.year}) should be the same."
    
    # store all output in dictionary
    vars_dict = {}

    # incomplete year of data in 1978
    if start_date.year == 1978:
        PPD_filename = 'icemotion_daily_nh_25km_19781101_19781231_v4.1.nc'

    # open file corresponding to date's year
    PP_mainpath = PPD_drift_path+PPD_filename
    PP_file = PP_mainpath.format(start_date.year,start_date.year)
    ds = xr.open_dataset(PP_file)
    ds.close()

    # cartopy projection for xx, yy, u_EASE, v_EASE
    vars_dict['proj'] = grab_projinfo_PPdrift(ds)

    
    # convert data times to CFTimeIndex
    CFTime = xr.cftime_range(f'{ds.time.values[0].year}', periods=len(ds.time), calendar="julian")
    # check that all dates match
    (CFTime.values).all() == (ds.time.values).all()
    # assign time variable to ds
    ds['time'] = CFTime.to_datetimeindex()
    
    # crop data to desired date range
    ds_date = ds.sel(time=slice(start_date, end_date))
    # convert to datetime
    ds_date['time'] = pd.to_datetime(ds_date.time.values)
        
        
    vars_dict['ds'] = ds_date

    # projected drift components
    u_EASE = ds_date.u.values
    v_EASE = ds_date.v.values

    # projected coords
    xx, yy = np.meshgrid(ds_date.x.values, ds_date.y.values)

    # lat/lon coords
    lat = ds_date.latitude.values
    lon = ds_date.longitude.values
    lon[lon<0]+=360

    # crop data to lat/long range
    lon, lat, [xx, yy, u_EASE, v_EASE] = crop_PPD_data(lon = lon, lat = lat, VARS = [xx, yy, u_EASE, v_EASE], lon_range=lon_range, lat_range=lat_range)

    # add to dict
    vars_dict['lon'] = lon
    vars_dict['lat'] = lat
    vars_dict['xx'] = xx
    vars_dict['yy'] = yy
    vars_dict['u_EASE'] = u_EASE
    vars_dict['v_EASE'] = v_EASE

    # if converted u, v desired, crop them
    if 'u' in return_vars or 'v' in return_vars:
        u, v = convert_PPD_vectors(lon = lon, u_EASE = u_EASE, v_EASE = v_EASE, fill_to_nan = True)
        vars_dict['u'] = u
        vars_dict['v'] = v

    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    if len(return_data) == 1:
        return_data = return_data[0]
    return return_data

#////////////////////////////
#  grab_projinfo_PPdrift ///
#//////////////////////////
#---------------------------------------------------------------------
# create cartopy projection for x, y coordinates of EASE grid for PPD data
#---------------------------------------------------------------------
# DEPENDENCIES:
import xarray as xr
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
#---------------------------------------------------------------------
def grab_projinfo_PPdrift(ds, quiet = True):
    
    """Grab projection info from NSIDC PP sea ice drift (NSIDC-0116) data (doi: 10.5067/MPYG15WAA4WX)

INPUT: 
- ds: sea ice drift data opened with xarray
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
- ice_projection: cartopy projection from data projection info

DEPENDENCIES:
import xarray as xr
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs

Latest recorded update:
06-28-2022
    """
    
    spatial = ds.crs
    
    # grab parameters from crs spatial attributes
    semimajor = float(spatial.proj4text[spatial.proj4text.find('+a=')+3:].split(' ')[0])
    semiminor = float(spatial.proj4text[spatial.proj4text.find('+b=')+3:].split(' ')[0])
    central_longitude = float(spatial.proj4text[spatial.proj4text.find('+lon_0=')+7:].split(' ')[0])
    central_latitude = float(spatial.proj4text[spatial.proj4text.find('+lat_0=')+7:].split(' ')[0])

    if not quiet:
        print(f'>>> data provided in {spatial.grid_mapping_name} projection from the {spatial.long_name}')
        print(f'  - semi_major_axis: {semimajor}')
        print(f'  - semi_minor_axis: {semiminor}')
        print(f'  - central_longitude: {central_longitude}')
        print(f'  - central_latitude: {central_latitude}')
        print(f'  - proj4text: {spatial.proj4text}')

    # create ice projection from info
    ice_projection = ccrs.LambertAzimuthalEqualArea(central_longitude=central_longitude, 
                                           central_latitude=central_latitude,
                                           globe=ccrs.Globe(semimajor_axis = semimajor, semiminor_axis = semiminor))
    
    return ice_projection