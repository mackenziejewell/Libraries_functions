# DEPENDENCIES:
import os
from datetime import datetime, timedelta
import xarray as xr
import netCDF4
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
from metpy.units import units
from pyhdf.SD  import *

def grab_projinfo_SIC(ds, quiet = True):

    """Grab projection info from Uni B AMSR2-AMSRE sea ice concentration data (doi: 10.1029/2005JC003384)
    https://seaice.uni-bremen.de/sea-ice-concentration/amsre-amsr2/

INPUT: 
- ds: data opened with xarray
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
- projection: cartopy projection from data projection info

DEPENDENCIES:
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
import xarray as xr

Latest recorded update:
01-26-2024
    """
    
    
    # grab projection info
    #---------------------
    CRS = ds.polar_stereographic.attrs

    # grab parameters from crs spatial attributes
    central_meridian = int(CRS['straight_vertical_longitude_from_pole'])
    semimajor = CRS['semi_major_axis']
    inv_flat = CRS['inverse_flattening']
    standard_parallel = int(CRS['standard_parallel'])
    
    if quiet != True:
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





def grab_ASI_SIC(date = datetime(year = 2015, month = 3, day = 24), 
                 main_path = '/Volumes/Jewell_EasyStore/UniB-ASI-SIC-n6250/',
                 coord_file = 'LongitudeLatitudeGrid-n6250-Arctic.hdf',
                 hemisphere = 'n',
                 resolution = '6250',
                 version = 'v5.4',
                 return_vars = ['xx', 'yy', 'lon', 'lat', 'sic', 'proj', 'ds'], 
                 include_units = False,
                 annual_folders = True, 
                 return_dict = False,
                 quiet = True):
    
    """Grab projection info from Uni B AMSR2-AMSRE sea ice concentration data (doi: 10.1029/2005JC003384)
    https://seaice.uni-bremen.de/sea-ice-concentration/amsre-amsr2/
    
    If annual_folders == True, expects the following directory stucture:
    
    main_path (top directory)
    │
    ├─── coord_file.hdf
    │
    ├─── sub_folder (by year, e.g. asi-AMSR2-n6250-2012)
    │    ├── file1.hdf (by day, e.g. asi-AMSR2-n6250-20120705-v5.4.nc)
    │    └── file2.hdf (by day)

    Else if annual_folders == False, expects all files in main folder:

    main_path (top directory)
    │
    ├─── coord_file.hdf
    ├─── file1.hdf (by day, e.g. asi-AMSR2-n6250-20120705-v5.4.nc)
    ├─── file2.hdf (by day)
    
    
    
INPUT: 
- date: date of data to open (datetime object)
    default: datetime(year = 2015, month = 3, day = 24)
- main_path: path to top directory where SIC data are stored (string)
    default: main_path = '/Volumes/Jewell_EasyStore/UniB-ASI-SIC/'
- coord_file: name of lat-lon coordinate file matching data resolution / hemisphere.
    Must be located in main_path!
    (default is for 6.25 km Arctic: 'LongitudeLatitudeGrid-n6250-Arctic.hdf')
    (file from https://data.seaice.uni-bremen.de/grid_coordinates/)
info used to determine file name and subfolder where data is stored:
- hemisphere: hemisphere of data as string, using UniBremen convention, either 'n' (north) or 's' (south)
- resolution: resolution of data as string, using UniBremen convention (default: '6250' for 6.25 km)
- version: version of data as string, using UniBremen convention (default: 'v5.4', most recent as of Jan 2024)
- return_vars: variables/attributes to return in specified order. (list)
    Can include any or all OUTPUT variables in any order: 
    default: ['xx', 'yy', 'lon', 'lat', 'sic', 'proj', 'ds']
- include_units: bool, whether or not to include units for relevant variables
- annual_folders: bool, whether or not files stored in annual subdirectories (default: True)
- return_dict: bool If True, return data as dictionary. If False, return list matching order of return_vars. (default: True)
- quiet: bool, whether or not to supress print statements (default: True)

OUTPUT:
List of any or all of variables specified in return_vars:
- xx: x-coordinates of projected SIC data (M x N array)
- yy: y-coordinates of projected SIC data (M x N array)
- lon: longitudes of projected data (M x N array)
- lat: latitudes of projected data (M x N array)
- sic: SIC data (M x N array)
- proj: cartopy projection from data projection info
- ds: xarray data frame containing opened SIC data

DEPENDENCIES:
import os
from datetime import datetime, timedelta
import xarray as xr
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
from metpy.units import units
from pyhdf.SD  import *
import netCDF4 (might be needed for xr to open ds?)

# homemade function:
grab_projinfo_SIC

Latest recorded update:
10-10-2024
    """

    # import numpy as np
    # import os
    # import xarray as xr
    # from datetime import datetime, timedelta
    # from metpy.units import units
    
    # assert input variable types
#     assert str(type(date)) == "<class 'datetime.datetime'>", f'date must be datetime object, not {type(date)}'
    assert type(return_vars) == list, f'return_vars must be list, not {type(return_vars)}'
    assert type(main_path) == str, f'main_path must be string, not {type(main_path)}'
    assert type(quiet) == bool, f'suppress_prints must be bool, not {type(suppress_prints)}'
    
    assert len(return_vars) > 0, 'return_vars list is empty. Must have length >=1'
    assert os.path.isdir(main_path), f'"{main_path}" not an existing directory'
    
    if not quiet:
        print('Provided input\n--------------')
        print(f'hemisphere: {hemisphere}')
        print(f'resolution: {resolution}')
        print(f'version: {version}')
        print(f'date: {date.date()}\n')
        
    # first check whether date is within data gap
    #--------------------------------------------
    # no data during gap from Oct 4 2011 - July 2 2012 
    # last and first dates with complete data outside gap
    gap_i_date = datetime(2011, 10, 3, 0, 0)
    gap_f_date = datetime(2012, 7, 3, 0, 0)

    # determine whether date is in data gap
    if (gap_i_date < date) & (date < gap_f_date):
        if not quiet:
            print(f"!!!\nWARNING >>> {date.date()} is within Oct 2011 - July 2012 data gap \n!!!\n")
    
    # create subfolder and file name
    #-------------------------------
    # grab date information
    datestring = date.strftime('%Y%m%d')
    year = date.year
    
    # determine which platform for subfolder and file
    if year <= 2011:
        platform = 'AMSRE'
    else:
        platform = 'AMSR2'
    
    # generate subfolder name 
    # e.g. 'asi-AMSR2-n6250-2012'
    if annual_folders:
        sub_folder = 'asi-{}-{}{}-{}/'.format(platform, hemisphere, resolution, year)
    else:
        sub_folder = ''
        
    assert os.path.isdir(main_path + sub_folder), f'"{main_path + sub_folder}" not an existing directory'
    if not quiet:
        print(f'main folder: {main_path}')
        print(f'sub folder: {sub_folder}')
    
    # generate file name 
    # e.g. 'asi-AMSR2-n6250-20120704-v5.4.nc' or 'asi-n6250-20120704-v5.4.nc' for AMSRE files
    if str(platform) == 'AMSRE':
        file_name = 'asi-{}{}-{}-{}.nc'.format(hemisphere, resolution, datestring, version)
    else:
        file_name = 'asi-{}-{}{}-{}-{}.nc'.format(platform, hemisphere, resolution, datestring, version)
    data_path = main_path + sub_folder + file_name
    assert os.path.isfile(data_path), f'"{file_name}" not an existing file in {sub_folder}'
    if not quiet:
        print(f'file name: {file_name}')

    # actually open data set
    ds = xr.open_dataset(data_path)
    ds.close()
    
    # store all output in dictionary
    vars_dict = {}
    
    # use grab_proj_info function to grab projection attributes and geo grid
    vars_dict['proj'] = grab_projinfo_SIC(ds)

    # netcdf4 dataset
    vars_dict['ds'] = ds
    
    
    if include_units:
        # sic
        vars_dict['sic'] = ds.z.values*units('percent')
        # projected geographic coordinates
        vars_dict['xx'], vars_dict['yy'] = np.meshgrid(ds.x.values*units(ds.x.units), ds.y.values*units(ds.y.units))
        
    else:
        vars_dict['sic'] = ds.z.values
        vars_dict['xx'], vars_dict['yy'] = np.meshgrid(ds.x.values, ds.y.values)

    
    if 'lon' in return_vars or 'lat' in return_vars:

        # grab lat / lon coordinates from coord_file
        assert os.path.isfile(main_path+coord_file), f'"{coord_file}" not an existing file in {main_path}'
        if not quiet:
            print(f'coordinate file: {coord_file}')
        
        # open file
        try:
            f = SD(main_path+coord_file, mode=1) 
        # raise an error if file can't be opened 
        except Exception as e:
            print(f"{e}, error opening {main_path+coord_file}")
        lon = f.select('Longitudes')[:]
        lat = f.select('Latitudes')[:]

        if include_units:
            # lats / lons
            vars_dict['lon'] = lon*units('degreeE')
            vars_dict['lat'] = lat*units('degreeN')
        else:
            vars_dict['lon'] = lon
            vars_dict['lat'] = lat
        

    # return as dictionary
    if return_dict:
        return vars_dict

    # return as list matching order of return_vars
    else:
        # save specified variables to list for output
        return_data = [vars_dict[var] for var in return_vars]
    
        if len(return_data) == 1:
            return_data = return_data[0]
    
        return return_data
        
