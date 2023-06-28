# DEPENDENCIES:
import os
from datetime import datetime
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc

def grab_projinfo_TB(ds, suppress_prints = True):
    
    """Grab projection info from NSIDC brightness temperature data (NSIDC-0001, doi: 10.5067/MXJL42WSXTS1)

INPUT: 
- ds: data opened with netCDF4
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
- projection: cartopy projection from data projection info

DEPENDENCIES:
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc

Latest recorded update:
04-25-2023
    """
    
    spatial = ds.variables['crs']
#     spatial = ds.crs
    
    # grab parameters from crs spatial attributes
    semimajor = spatial.semi_major_axis
    semiminor = float(spatial.proj4text[spatial.proj4text.find('+b=')+3:].split(' ')[0])
    inv_flat = spatial.inverse_flattening
    central_meridian = int(spatial.straight_vertical_longitude_from_pole)
    latitude_of_origin = int(spatial.standard_parallel)

    if suppress_prints != True:
        print(f'>>> data provided in {spatial.grid_mapping_name} projection from the {spatial.long_name}')
        print(f'  - semi_major_axis: {semimajor}')
        print(f'  - semi_minor_axis: {semiminor}')
        print(f'  - inverse_flattening: {inv_flat}')
        print(f'  - straight_vertical_longitude_from_pole: {central_meridian}')
        print(f'  - standard_parallel: {latitude_of_origin}')
        print(f'  - proj4text: {spatial.proj4text}')

    # create projection from info
    projection = ccrs.NorthPolarStereo(central_longitude=int(central_meridian), 
                                           globe=ccrs.Globe(semimajor_axis = semimajor, semiminor_axis = semiminor, 
                                                            inverse_flattening=inv_flat),
                                           true_scale_latitude=int(latitude_of_origin))
    
    return projection



def grab_TB_NSIDC0001(date = datetime(year = 2015, month = 3, day = 24), 
                      Tb_datapath= '/Volumes/Jewell_EasyStore/NSIDC-0001_TB/', 
                      Tb_name = 'NSIDC0001_TB_PS_N25km_{}_v6.0.nc',
                      return_vars = ['xx', 'yy', 'TB_19H', 'TB_19V', 'TB_37H', 'TB_37V', 'proj', 'ds'], 
                      throw_error_miss = True,
                      suppress_prints = True):
   
    
    """Grab daily brightness temperatures (Tb) from NSIDC data (NSIDC-0001, doi: 10.5067/MXJL42WSXTS1)

INPUT: 
- date: date of data to open (datetime object)
    default: datetime(year = 2015, month = 3, day = 24)
- Tb_datapath: path to directory where SIC data are stored (string)
    default: Tb_datapath = '/Volumes/Jewell_EasyStore/Tb/'
- Tb_name = naming convention of sic data files with {} indicating location of date (YYYYmmdd) in string (string)
    default: 'NSIDC0001_TB_PS_N25km_{}_v6.0.nc'
- return_vars: variables/attributes to return in specified order. (list)
    Can include any or all OUTPUT variables in any order: 
    default: ['xx', 'yy', 'TB_19H', 'TB_19V', 'TB_37H', 'TB_37V', 'PR_19', 'PR_37', 'GR_37_19_V', 'GR_37_19_H', 'proj', 'ds']
- throw_error_miss: bool, whether or not to throw an error if data not in dataset, e.g. 2/20/2021 (default: True)
    If False, return empty string for all vars besides 'xx', 'yy', 'proj', and 'ds' in dictionary when data missing
    If True, raise exception. 
- suppress_prints: bool, whether or not to supress print statements (default: True)

OUTPUT:
List of any or all of variables specified in return_vars:
- xx: x-coordinates of projected SIC data (M x N array)
- yy: y-coordinates of projected SIC data (M x N array)
- TB_19H: Horizontal polarization of 19 Hz Tb data (M x N array)
- TB_19V: Vertical polarization of 19 Hz Tb data (M x N array)
- TB_37H: Horizontal polarization of 37 Hz Tb data (M x N array)
- TB_37V: Vertical polarization of 37 Hz Tb data (M x N array)
- PR_19: polarization ratio @ 19Ghz as defined in Cavalieri et al., (1984) doi:10.1029/JD089iD04p05355
- PR_37: polarization ratio @ 37Ghz as defined in Cavalieri et al., (1984) doi:10.1029/JD089iD04p05355
- GR_37_19_V: gradient ratio w/ V polarization as defined in Cavalieri et al., (1984) doi:10.1029/JD089iD04p05355
- GR_37_19_H: gradient ratio w/ H polarization as defined in Cavalieri et al., (1984) doi:10.1029/JD089iD04p05355
- proj: cartopy projection from data projection info
- ds: xarray data frame containing opened Tb data (does not actually show Tb data, must be opened with NetCDF4)

DEPENDENCIES:
import os
from datetime import datetime
import xarray as xr
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc

# homemade function:
grab_projinfo_TB

Latest recorded update:
06-27-2023
    """
    
    # assert input variable types
    assert str(type(date)) == "<class 'datetime.datetime'>", f'date must be datetime object, not {type(date)}'
    assert type(return_vars) == list, f'return_vars must be list, not {type(return_vars)}'
    assert type(Tb_name) == str, f'Tb_name must be string, not {type(Tb_name)}'
    assert type(Tb_datapath) == str, f'Tb_datapath must be string, not {type(Tb_datapath)}'
    assert type(suppress_prints) == bool, f'suppress_prints must be bool, not {type(suppress_prints)}'
    
    # open Tb data
    filename = Tb_name.format(date.strftime('%Y%m%d'))
    data_path = Tb_datapath+filename
    
    if suppress_prints == False:
        print(f' >>> opening {data_path}')
    
    assert len(return_vars) > 0, 'return_vars list is empty. Must have length >=1'
    assert os.path.isdir(Tb_datapath), f'"{Tb_datapath}" not an existing directory'
    assert os.path.isfile(data_path), f'"{filename}" not an existing file in {Tb_datapath}'
    
    
    # open data with xarray to grab time / projection info
#     ds = xr.open_dataset(data_path)
#     ds.close()
    
    # remove time dimension since data are stored in single-time slices
#     ds_time = ds.time[0]
#     if suppress_prints == False:
#         print(f' time: {ds_time.values}')
#     ds = ds.sel(time = ds_time)

    
    # open data with netCDF4 to import actual data
    ds = nc.Dataset(data_path)

    # store all output in dictionary
    vars_dict = {}
    
    # use grab_proj_info function to grab projection attributes and geo grid
    vars_dict['proj'] = grab_projinfo_TB(ds)
    
    # projected geographic coordinates
    vars_dict['xx'], vars_dict['yy'] = np.meshgrid(ds.variables['x'][:], ds.variables['y'][:])

    # netcdf4 dataset
    vars_dict['ds'] = ds
    
    
    # for var in ds.variables.values():
    #     print(var)
    # for var in ds.groups.values():
    #     print(var)
#     for var in ds.groups['F17'].variables:
#             print(var)    
    
    # if sensor data exists, default to first sensor data (for example F17, if both F17, and F18 data are stored)
    # could later add option for user to set preference for sensor platform if available on given date
    if len(ds.groups) > 0:
 
        # first sensor name
        # e.g. if sat == 'F17', this would grab ds.groups[F17']['TB_F17_19H']
        sat = list(ds.groups)[0]

        if suppress_prints == False:
            print(f'satellites: {list(ds.groups)}')
            print(f'select {sat}')

        vars_dict['TB_19H'] = ds.groups[sat][f'TB_{sat}_19H'][0,:,:]
        vars_dict['TB_19V'] = ds.groups[sat][f'TB_{sat}_19V'][0,:,:]
        vars_dict['TB_37H'] = ds.groups[sat][f'TB_{sat}_37H'][0,:,:]
        vars_dict['TB_37V'] = ds.groups[sat][f'TB_{sat}_37V'][0,:,:]   

        # polarization and gradient ratios as in Cavalieri et al. (1984)
        vars_dict['PR_19'] = (vars_dict['TB_19V'] - vars_dict['TB_19H']) / (vars_dict['TB_19V'] + vars_dict['TB_19H'])
        vars_dict['PR_37'] = (vars_dict['TB_37V'] - vars_dict['TB_37H']) / (vars_dict['TB_37V'] + vars_dict['TB_37H'])
        vars_dict['GR_37_19_V'] = (vars_dict['TB_37V'] - vars_dict['TB_19V']) / (vars_dict['TB_37V'] + vars_dict['TB_19V'])
        vars_dict['GR_37_19_H'] = (vars_dict['TB_37H'] - vars_dict['TB_19H']) / (vars_dict['TB_37H'] + vars_dict['TB_19H'])
    
    # if sensor data does not exist, either throw error or assign as empty
    else:
        if throw_error_miss:
            raise Exception(f"no sensor data variable found in ds groups for {date.strftime('%b %d, %Y')}")
        else:
            for VAR in ['TB_19H', 'TB_19V', 'TB_37H', 'TB_37V', 'PR_19', 'PR_37', 'GR_37_19_V', 'GR_37_19_H']:
                vars_dict[VAR] = ''

    
    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    
    if len(return_data) == 1:
        return_data = return_data[0]

    return return_data
