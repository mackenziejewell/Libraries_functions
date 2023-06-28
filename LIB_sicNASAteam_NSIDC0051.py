# DEPENDENCIES:
import os
from datetime import datetime
import xarray as xr
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs


def grab_projinfo_sicNASAteam(ds, quiet = True):
    
    """Grab projection info from NASA team (NSIDC-0051) sea ice concentration data (doi: 10.5067/MPYG15WAA4WX)

INPUT: 
- ds: sea ice concentration data opened with xarray
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
- ice_projection: cartopy projection from sea ice concentration data projection info

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
    semimajor = spatial.semi_major_axis
    semiminor = float(spatial.proj4text[spatial.proj4text.find('+b=')+3:].split(' ')[0])
#     inv_flat = spatial.inverse_flattening
    central_meridian = int(spatial.straight_vertical_longitude_from_pole)
    latitude_of_origin = int(spatial.standard_parallel)

    if not quiet:
        print(f'>>> data provided in {spatial.grid_mapping_name} projection from the {spatial.long_name}')
        print(f'  - semi_major_axis: {semimajor}')
        print(f'  - semi_minor_axis: {semiminor}')
        print(f'  - inverse_flattening: {inv_flat}')
        print(f'  - straight_vertical_longitude_from_pole: {central_meridian}')
        print(f'  - standard_parallel: {latitude_of_origin}')
        print(f'  - proj4text: {spatial.proj4text}')

    # create ice projection from info
    ice_projection = ccrs.NorthPolarStereo(central_longitude=int(central_meridian), 
                                           globe=ccrs.Globe(semimajor_axis = semimajor, semiminor_axis = semiminor),
#                                                             ,inverse_flattening=inv_flat),
                                           true_scale_latitude=int(latitude_of_origin))
    
    return ice_projection




def grab_sicNASAteam(date = datetime(year = 2000, month = 1, day = 1), 
                     sic_datapath= '/Volumes/Jewell_EasyStore/NSIDC-0051_NTsic/', 
                     SIC_name = 'NSIDC0051_SEAICE_PS_N25km_{}_v2.0.nc',
                     return_vars = ['xx', 'yy', 'sic', 'proj', 'ds'], 
                     throw_error_miss = True, quiet = True):
    
    """Grab daily sea ice concentrations from NASA team (NSIDC-0051) sea ice concentration data (doi: 10.5067/MPYG15WAA4WX)

INPUT: 
- date: date of SIC data to open (datetime object)
    default: datetime(year = 2015, month = 3, day = 24)
- sic_datapath: path to directory where SIC data are stored (string)
    default: sic_datapath = '/Volumes/Jewell_EasyStore/NSIDC-0051_NTsic/'
- SIC_name = naming convention of sic data files with {} indicating location of date (YYYYmmdd) in string (string)
    default: 'NSIDC0051_SEAICE_PS_N25km_{}_v2.0.nc'
- return_vars: variables/attributes to return in specified order. (list)
    Can include any or all OUTPUT variables in any order: 
    default: ['xx', 'yy', 'sic', 'proj', 'ds']
- throw_error_miss: bool, whether or not to throw an error if SIC data not in dataset, e.g. late 1987 (default: True)
    If False, return empty string for 'sic' in dictionary when data missing
    If True, raise exception. 
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT: 
- xx: x-coordinates of projected SIC data (M x N array)
- yy: y-coordinates of projected SIC data (M x N array)
- sic: SIC data (M x N array)
- proj: cartopy projection from sea ice concentration data projection info
- ds: xarray data frame containing opened SIC data

DEPENDENCIES:
import os
from datetime import datetime
import xarray as xr
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs

# homemade function:
grab_projinfo_sicNASAteam

Latest recorded update:
06-28-2022
    """
    
    # assert input variable types
    assert str(type(date)) == "<class 'datetime.datetime'>", f'date must be datetime object, not {type(date)}'
    assert type(return_vars) == list, f'return_vars must be list, not {type(return_vars)}'
    assert type(SIC_name) == str, f'SIC_name must be string, not {type(SIC_name)}'
    assert type(sic_datapath) == str, f'sic_datapath must be string, not {type(sic_datapath)}'
    assert type(quiet) == bool, f'quiet must be bool, not {type(quiet)}'
    
    # open SIC data
    filename = SIC_name.format(date.strftime('%Y%m%d'))
    data_path = sic_datapath+filename
    
    if not quiet:
        print(f' >>> opening {data_path}')
    
    assert len(return_vars) > 0, 'return_vars list is empty. Must have length >=1'
    assert os.path.isdir(sic_datapath), f'"{sic_datapath}" not an existing directory'
    assert os.path.isfile(data_path), f'"{filename}" not an existing file in {sic_datapath}'
    
    
    # open data
    ds = xr.open_dataset(data_path)
    ds.close()
    
    # remove time dimension since data are stored in single-time slices
    ds_time = ds.time[0]
    if not quiet:
        print(f' time: {ds_time.values}')
    ds = ds.sel(time = ds_time)

    # store all output in dictionary
    vars_dict = {}
    
    vars_dict['ds'] = ds
    
    # use grab_proj_info function to grab projection attributes and geo grid
    vars_dict['proj'] = grab_projinfo_sicNASAteam(ds)
    vars_dict['xx'], vars_dict['yy'] = np.meshgrid(ds.x.values, ds.y.values)
    
    # Grab SIC
    #---------------------------------------------------------------
    # should be 5 variables when SIC data avilable
    # crs, x, y, time, and the name of the SIC data variable
    # 'F17_ICECON' for 2008 - present
    # 'F13_ICECON' for 1995 - 2007
    # 'F11_ICECON' for 1991 - 1995
    # 'F08_ICECON' for 1987 - 1991 (none for late 1987, data gap)
    # 'N07_ICECON' 1978 - 1987
    
    # list all keys, which include coords and SIC variable name if exists
    all_keys = list(ds.variables)
    coords = ['crs', 'x', 'y', 'time']

    # look through all_keys to find name of SIC variable. 
    # If exists, replace sic_variable. If not, sic_variable empty
    sic_variable = ''
    for key in all_keys:
        if key not in coords:
            sic_variable = key

    # if sic_variable exists in ds, grab SIC
    if len(sic_variable) > 0:
        vars_dict['sic'] = ma.masked_where(ds.get(sic_variable)>1, ds.get(sic_variable))
    # if sic_variable does not exist, either throw error or assign as None
    else:
        if throw_error_miss:
            raise Exception(f"no SIC data variable found in ds variables: {all_keys}")
        else:
            vars_dict['sic'] = ''
    #---------------------------------------------------------------

    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    
    if len(return_data) == 1:
        return_data = return_data[0]

    return return_data
