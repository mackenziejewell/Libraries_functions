# DEPENDENCIES:
import os
from datetime import datetime
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
import xarray as xr

def grab_projinfo_SIC(ds, suppress_prints = True):

    """Grab projection info from 1 km MODIS-AMSR2 sea ice concentration data (doi: 10.3390/rs12193183)
    https://seaice.uni-bremen.de/sea-ice-concentration/modis-amsr2/

INPUT: 
- ds: data opened with xarray
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- projection: cartopy projection from data projection info

DEPENDENCIES:
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs

Latest recorded update:
10-01-2023
    """
    
    
    # grab projection info
    #---------------------
    CRS = ds.polar_stereographic.attrs

    # grab parameters from crs spatial attributes
    central_meridian = int(CRS['straight_vertical_longitude_from_pole'])
    semimajor = CRS['semi_major_axis']
    inv_flat = CRS['inverse_flattening']
    standard_parallel = int(CRS['standard_parallel'])
    
    if suppress_prints != True:
        print(f'>>> data provided in polar_stereographic projection')
        print(f'  - semi_major_axis: {semimajor}')
        print(f'  - inverse_flattening: {inv_flat}')
        print(f'  - straight_vertical_longitude_from_pole: {central_meridian}')
        print(f'  - standard_parallel: {standard_parallel}')
        print(f'  - proj4text: {ds.attrs["proj4string"]}')
        
    # create projection from info
    projection = ccrs.NorthPolarStereo(central_longitude=central_meridian, 
                                           globe=ccrs.Globe(semimajor_axis = semimajor,inverse_flattening=inv_flat),
                                           true_scale_latitude=standard_parallel)

    return projection


def grab_SIC_MODISAMSR(date = datetime(year = 2015, month = 3, day = 24), 
                      file_datapath= '/Volumes/Seagate_Jewell/KenzieStuff/SIC_MODIS_AMSR2/', 
                      file_name = 'sic_modis-aqua_amsr2-gcom-w1_merged_nh_1000m_{}.nc',
                      geo_file = 'coordinates_npstere_1km_arctic.nc',
                      return_vars = ['xx', 'yy', 'lon', 'lat', 'sic_merged', 'unc_sic_merged', 'sic_modis', 'sic_amsr2', 'proj', 'ds'], 
                      suppress_prints = True):
   
    
    """Grab daily sea ice concentration (SIC) from 1 km MODIS-AMSR2 SIC data (doi: 10.3390/rs12193183)
    https://seaice.uni-bremen.de/sea-ice-concentration/modis-amsr2/

INPUT: 
- date: date of data to open (datetime object)
    default: datetime(year = 2015, month = 3, day = 24)
- file_datapath: path to directory where SIC data are stored (string)
    default: file_datapath = '/Volumes/Seagate_Jewell/KenzieStuff/SIC_MODIS_AMSR2/'
- file_name = naming convention of sic data files with {} indicating location of date as YYYYmmdd
    default: 'sic_modis-aqua_amsr2-gcom-w1_merged_nh_1000m_{}.nc'
- geo_file: name of file listing latitude/longitude coordinates matching xx, yy grid of projected SIC file.
    default: 'coordinates_npstere_1km_arctic.nc'
    (file from https://data.seaice.uni-bremen.de/modis_amsr2/netcdf/)
- return_vars: variables/attributes to return in specified order. (list)
    Can include any or all OUTPUT variables in any order: 
    default: ['xx', 'yy', 'lon', 'lat', 'sic_merged', 'unc_sic_merged', 'sic_modis', 'sic_amsr2', 'proj', 'ds']
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
List of any or all of variables specified in return_vars:
- xx: x-coordinates of projected SIC data (M x N array)
- yy: y-coordinates of projected SIC data (M x N array)
- lon: longitudes of projected data (M x N array)
- lat: latitudes of projected data (M x N array)
- sic_merged: merged MODIS-AMSR2 SIC data (M x N array)
- unc_sic_merged: uncertainty in merged MODIS-AMSR2 SIC data (M x N array)
- sic_modis: SIC estimate from MODIS data only (M x N array)
- sic_amsr2: SIC estimate from AMSR2 data only (M x N array)
- proj: cartopy projection from data projection info
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
grab_projinfo_SIC

Latest recorded update:
10-01-2023
    """
    
    # assert input variable types
    assert str(type(date)) == "<class 'datetime.datetime'>", f'date must be datetime object, not {type(date)}'
    assert type(return_vars) == list, f'return_vars must be list, not {type(return_vars)}'
    assert type(file_name) == str, f'file_name must be string, not {type(file_name)}'
    assert type(file_datapath) == str, f'file_datapath must be string, not {type(file_datapath)}'
    assert type(suppress_prints) == bool, f'suppress_prints must be bool, not {type(suppress_prints)}'
    
    # open SIC data
    filename = file_name.format(date.strftime('%Y%m%d'))
    data_path = file_datapath+filename
    
    if not suppress_prints:
        print(f' >>> opening {data_path}')
    
    assert len(return_vars) > 0, 'return_vars list is empty. Must have length >=1'
    assert os.path.isdir(file_datapath), f'"{file_datapath}" not an existing directory'
    assert os.path.isfile(data_path), f'"{filename}" not an existing file in {file_datapath}'
    
    
    # open data with xarray to grab projection info
    ds = xr.open_dataset(data_path)
    ds.close()

    # store all output in dictionary
    vars_dict = {}
    
    # use grab_proj_info function to grab projection attributes and geo grid
    vars_dict['proj'] = grab_projinfo_SIC(ds)
   
    # projected geographic coordinates
    vars_dict['xx'], vars_dict['yy'] = np.meshgrid(ds.x.values, ds.y.values)

    # netcdf4 dataset
    vars_dict['ds'] = ds
    
    # grab SIC data 
    # (for some reason dimension 1 must be reversed to match x/y and lat/lon)
    if 'sic_merged' in return_vars:
        vars_dict['sic_merged'] = ds.sic_merged.values[::-1,:]
    if 'unc_sic_merged' in return_vars:
        vars_dict['unc_sic_merged'] = ds.unc_sic_merged.values[::-1,:]
    if 'sic_modis' in return_vars:
        vars_dict['sic_modis'] = ds.sic_modis.values[::-1,:]
    if 'sic_amsr2' in return_vars:
        vars_dict['sic_amsr2'] = ds.sic_amsr2.values[::-1,:]
    
    
    # grab geo coordinates
    #---------------------
    if 'lat' in return_vars or 'lon' in return_vars:
        ds = xr.open_dataset(file_datapath+geo_file)
        ds.close()
        lon = ds.lon.values
        lon[lon<0]+=360
        lat = ds.lat.values

        vars_dict['lon'] = lon
        vars_dict['lat'] = lat
        

    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    
    if len(return_data) == 1:
        return_data = return_data[0]

    return return_data
