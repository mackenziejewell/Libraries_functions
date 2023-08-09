#/////////////////////////////
#  create_AMSRnc_cellarea ///
#///////////////////////////
#---------------------------------------------------------------------
# Create cropped pixel area file corresponding to locally-created .nc files
#---------------------------------------------------------------------
#/////////////////////////////
#  create_AMSRnc_landmask ///
#///////////////////////////
#---------------------------------------------------------------------
# Create cropped landmask corresponding to locally-created .nc files
#---------------------------------------------------------------------
#/////////////////////////////
#  convert_AMSR_he5_to_nc ///
#///////////////////////////
#---------------------------------------------------------------------
# Crop and convert locally-stored .he5 files to smaller .nc files
#---------------------------------------------------------------------
#////////////////////
#  read_AMSR_he5 ///
#//////////////////
#---------------------------------------------------------------------
# Grab 89 GHz brightness temperatures from locally-stored .he5 data
#---------------------------------------------------------------------
#////////////////
#  crop_data ///
#//////////////
#---------------------------------------------------------------------
# Crop AMSR data and coords to within given lat/lon range
#---------------------------------------------------------------------


#/////////////////////////////
#  create_AMSRnc_cellarea ///
#///////////////////////////
#---------------------------------------------------------------------
# Create cropped pixel area file corresponding to locally-created .nc files
#---------------------------------------------------------------------
# DEPENDENCIES:
from pyhdf.SD import SD, SDC
import datetime
import numpy as np
import xarray as xr
# from LIB_AMSR89_NSIDCAU_SI6 import read_AMSR_he5, crop_data
#---------------------------------------------------------------------
def create_AMSRnc_cellarea(AMSR_area_file = '/Users/mackenziejewell/Data/AMSR89Tb_NSIDCAU_SI6/NSIDC0771_CellArea_PS_N6.25km_v1.0.nc',
                           AMSR_he5_file = '/Users/mackenziejewell/Data/AMSR89Tb_NSIDCAU_SI6/AMSR_U2_L3_SeaIce6km_B04_20210101.he5',
                           lon_range = [180, 260],
                           lat_range = [65, 90]):
    """Create cropped landmask corresponding to locally-created cropped 6.25 km 89 GHz Tb from locally-stored .he5 NSIDC data (NSIDC-AU_SI6, doi: 10.5067/NX1R09ORNOZN)
INPUT: 
- AMSR_area_file: path to locally stored NSIDC AMSR 6.25 km pixel area file (NSIDC0771_CellArea_PS_N6.25km_v1.0.nc)
    Downloaded from: https://nsidc.org/data/nsidc-0771/versions/1
    [Polar Stereographic Ancillary Grid Information, Version 1 NSIDC-0771 doi:10.5067/N6INPBT8Y104]
- AMSR_he5_file = path to locally stored NSIDC AMSR 6.25 km he5 file (doi: 10.5067/NX1R09ORNOZN)
    for spatial reference of desired lat/long grid to crop to
    (default: '/Users/mackenziejewell/Data/AMSR89Tb_NSIDCAU_SI6/AMSR_U2_L3_SeaIce6km_B04_20210101.he5')
- lon_range: [lonmin, lonmax] to crop, values in range (0,360). default: [180, 260]
- lat_range: [latmin, latmax] to crop. default: [65,90]

OUTPUT:
- lon_crop: (M x N) cropped longitude grid
- lat_crop: (M x N) cropped latitude grid
- area_crop: (M x N) area (m^2) cropped to provided lon/lat range.

DEPENDENCIES:
from pyhdf.SD import SD, SDC
import datetime
import numpy as np
import xarray as xr

# homemade function
from LIB_AMSR89_NSIDCAU_SI6 import read_AMSR_he5, crop_data


Latest recorded update:
08-08-2023
    """
    
    # open pixel area data (1792 x 1216)
    #---------------------------------
    ds = xr.open_dataset(AMSR_area_file)
    area = ds.cell_area.values # meters ^ 2

    # open original.hef data (1792 x 1216)
    #-------------------------------------
    # extract file info to input into read_AMSR_he5 
    AMSR_filepath = AMSR_he5_file[:AMSR_he5_file.find(AMSR_he5_file.split('/')[-1])]
    AMSR_datestring = AMSR_he5_file.split('/')[-1].split('.')[0].split('_')[-1]
    AMSR_date = datetime.datetime.strptime(AMSR_datestring, '%Y%m%d')

    # grab lat lon data from .he5 file
    [LAT, LON] = read_AMSR_he5(date=AMSR_date, 
                               filepath=AMSR_filepath, 
                               filename_convention='AMSR_U2_L3_SeaIce6km_B04_{}{}{}.he5', 
                               hemisphere='N', add_units=False, 
                               return_vars=['lat', 'lon'], quiet=False)
    
    # crop pixel area to
    # match files saved to .nc cropped as (M x N)
    #-------------------------------------
    lon_crop, lat_crop, [area_crop] = crop_data(lon=LON, lat=LAT, VARS=[area], 
                                                     lon_range=lon_range, lat_range=lat_range)

    return lon_crop, lat_crop, area_crop



# # create netCDF data set
# #-----------------------
# ds = xr.Dataset(data_vars=dict(
#     pixelarea = (["x", "y"], area_crop, 
#                {"description":"pixel area", "units": 'meters^2'})),
# coords=dict(
#     longitude = (["x", "y"], lon_crop, {"units": 'degrees_east'}),
#     latitude  = (["x", "y"], lat_crop, {"units": 'degrees_north'}),
# ),
# attrs=dict(
#     description = (f"Pixel area corresponding to 6.25 km 89 GHz brightness temperatures, cropped ([{lat_range[0]},{lat_range[1]}], [{lon_range[0]},{lon_range[1]}]). (NSIDC-0771, doi: 10.5067/N6INPBT8Y104)"), 
#     ),
# )
# ds.to_netcdf(path='/Users/mackenziejewell/Data/AMSR89Tb_NSIDCAU_SI6/AMSR_6km_area.nc', mode='w', format="NETCDF4")


#/////////////////////////////
#  create_AMSRnc_landmask ///
#///////////////////////////
#---------------------------------------------------------------------
# Create cropped landmask corresponding to locally-created .nc files
#---------------------------------------------------------------------
# DEPENDENCIES:
from pyhdf.SD import SD, SDC
import datetime
import numpy as np
# from LIB_AMSR89_NSIDCAU_SI6 import read_AMSR_he5, crop_data
#---------------------------------------------------------------------
def create_AMSRnc_landmask(AMSR_mask_file = '/Users/mackenziejewell/Data/AMSR89Tb_NSIDCAU_SI6/amsr_gsfc_6n.hdf',
                           AMSR_he5_file = '/Users/mackenziejewell/Data/AMSR89Tb_NSIDCAU_SI6/AMSR_U2_L3_SeaIce6km_B04_20210101.he5',
                           lon_range = [180, 260],
                           lat_range = [65, 90]):
    """Create cropped landmask corresponding to locally-created cropped 6.25 km 89 GHz Tb from locally-stored .he5 NSIDC data (NSIDC-AU_SI6, doi: 10.5067/NX1R09ORNOZN)
INPUT: 
- AMSR_mask_file: path to locally stored NSIDC AMSR 6.25 km landmask file
    Downloaded from: https://nsidc.org/data/user-resources/help-center/does-nsidc-have-tools-extract-and-geolocate-polar-stereographic-data
- AMSR_he5_file = path to locally stored NSIDC AMSR 6.25 km he5 file (doi: 10.5067/NX1R09ORNOZN)
    for spatial reference of desired lat/long grid to crop to
    (default: '/Users/mackenziejewell/Data/AMSR89Tb_NSIDCAU_SI6/AMSR_U2_L3_SeaIce6km_B04_20210101.he5')
- lon_range: [lonmin, lonmax] to crop, values in range (0,360). default: [180, 260]
- lat_range: [latmin, latmax] to crop. default: [65,90]

OUTPUT:
- lon_crop: (M x N) cropped longitude grid
- lat_crop: (M x N) cropped latitude grid
- land_mask: (M x N) land mask cropped to provided lon/lat range.
    1 = land or coast. 0 = water.

DEPENDENCIES:
from pyhdf.SD import SD, SDC
import datetime
import numpy as np

# homemade function
from LIB_AMSR89_NSIDCAU_SI6 import read_AMSR_he5, crop_data


Latest recorded update:
08-08-2023
    """
    
    # open landmask data (1792 x 1216)
    #---------------------------------
    f = SD(AMSR_mask_file,SDC.READ)
    land_ma = f.select('landmask')[:]
    
    
    # Values are 0 (water), 3 (coast), 216 (land)
    # https://nsidc.org/data/user-resources/help-center/does-nsidc-have-tools-extract-and-geolocate-polar-stereographic-data

    # open original.hef data (1792 x 1216)
    #-------------------------------------
    # extract file info to input into read_AMSR_he5 
    AMSR_filepath = AMSR_he5_file[:AMSR_he5_file.find(AMSR_he5_file.split('/')[-1])]
    AMSR_datestring = AMSR_he5_file.split('/')[-1].split('.')[0].split('_')[-1]
    AMSR_date = datetime.datetime.strptime(AMSR_datestring, '%Y%m%d')

    # grab lat lon data from .he5 file
    [LAT, LON] = read_AMSR_he5(date=AMSR_date, 
                               filepath=AMSR_filepath, 
                               filename_convention='AMSR_U2_L3_SeaIce6km_B04_{}{}{}.he5', 
                               hemisphere='N', add_units=False, 
                               return_vars=['lat', 'lon'], quiet=False)
    
    # crop landmask as data 
    # match files saved to .nc cropped as (M x N)
    #-------------------------------------
    lon_crop, lat_crop, [land_mask_crop] = crop_data(lon=LON, lat=LAT, VARS=[land_ma], 
                                                     lon_range=lon_range, lat_range=lat_range)
    
    # turn to integers
    land_mask = land_mask_crop.astype(int)
    
    # set land (216) and coast (3) to value of 1
    land_mask[land_mask==216]=1
    land_mask[land_mask==3]=1
    
    return lon_crop, lat_crop, land_mask


# # create netCDF data set
# #-----------------------
# ds = xr.Dataset(data_vars=dict(
#     landmask = (["x", "y"], LM, 
#                {"description":"0: ocean, 1: land/coast"})),
# coords=dict(
#     longitude = (["x", "y"], loc, {"units": 'degrees_east'}),
#     latitude  = (["x", "y"], lac, {"units": 'degrees_north'}),
# ),
# attrs=dict(
#     description = (f"Land mask corresponding to 6.25 km 89 GHz brightness temperatures, cropped ([{lat_range[0]},{lat_range[1]}], [{lon_range[0]},{lon_range[1]}]). (NSIDC-AU_SI6, doi: 10.5067/NX1R09ORNOZN)"), 
#     ),
# )
# ds.to_netcdf(path='/Users/mackenziejewell/Data/AMSR89Tb_NSIDCAU_SI6/AMSR_6km_landmask.nc', mode='w', format="NETCDF4")



#/////////////////////////////
#  convert_AMSR_he5_to_nc ///
#///////////////////////////
#---------------------------------------------------------------------
# Crop and convert locally-stored .he5 files to smaller .nc files
#---------------------------------------------------------------------
# DEPENDENCIES:
import os
import h5py
import numpy as np
import xarray as xr
import datetime
#---------------------------------------------------------------------
def convert_AMSR_he5_to_nc(file = '/Users/mackenziejewell/Downloads/AMSR_U2_L3_SeaIce6km_B04_20000101.he5', 
                           hemisphere = 'N', lon_range=[180, 260], lat_range=[65,90], 
                           new_filename = 'AMSR_89Ghz_{}.nc', save_folder = 'converted',
                           move_he5_file = False, allow_overwrites = False,
                           return_new_filename = False,
                           quiet = True):
    """Grab 6.25 km 89 GHz brightness temperatures from locally-stored AMSR-2 NSIDC data (NSIDC-AU_SI6, doi: 10.5067/NX1R09ORNOZN)
    and crop then save as much smaller .nc file.
    
INPUT: 
- file: name of stored .he5 data file, including path to file
    (default: '/Users/mackenziejewell/Downloads/AMSR_U2_L3_SeaIce6km_B04_20000101.he5')
- hemisphere: string specifying hemispheric data to grab.
    either 'N' for north hemisphere data or 'S' for south hemisphere data
- lon_range: [lonmin, lonmax] to crop, values in range (0,360). default: [180, 260]
- lat_range: [latmin, latmax] to crop. default: [65,90]

- new_filename: name with which to save new file (with {} inserted where date info will be placed)
        (default: 'AMSR_89Ghz_{}.nc')
- save_folder: name of folder within file's path to save file to. Will then copy used .he5 file to same folder.
                Do not add any slashes to folder name. (default: 'converted').
                If folder does not yet exist, will create folder.
- move_he5_file: bool, whether or not to move .he5 file to save_folder alongside new .nc file
    (default: False)              
- allow_overwrites: bool, whether or not to allow file overwrites (default: False)
- return_new_filename: bool, whether or not to return name of new .nc file. 
    if True, return filename. 
    if False (default), don't return anything.

quiet = False



- quiet: bool, whether or not to hide print statements (default: True)
    

OUTPUT:
List of any or all of variables specified in return_vars:
(or list of empty lists if data file does not exist)
- lat: (MxN) grid of latitude data
- lon: (MxN) grid of longitude data (0-360)

- 89H_DAY: (MxN) grid of H pol 89 GHz daily average Tb
- 89V_DAY: (MxN) grid of V pol 89 GHz daily average Tb
- 89H_ASC: (MxN) grid of H pol 89 GHz daily average ascending Tb
- 89V_ASC: (MxN) grid of V pol 89 GHz daily average ascending Tb
- 89H_DSC: (MxN) grid of H pol 89 GHz daily average descending Tb
- 89V_DSC: (MxN) grid of V pol 89 GHz daily average descending Tb


DEPENDENCIES:
import os
import h5py
import numpy as np
import xarray as xr
import datetime

# homemade function:
crop_data

Latest recorded update:
08-02-2023
    """
    
    
    # create function which determines whether or not to write files
    #---------------------------------------------------------------
    def check_writefile(file, allow_overwrites = allow_overwrites, quiet = quiet):

        # return True if file doesn't exist or it does but overwrites are allowed
        # return False is file exists and overwrites not allowed

        # check whether the file exists
        #------------------------------
        isFile = os.path.isfile(file)

        # if file already exists, check if it should be overwritten
        if isFile:
            if not quiet:
                print(f'\n!!! file already exists: {file} !!!')

            # determine whether or not to write file
            if allow_overwrites:
                if not quiet:
                    print('>>> Overwriting existing file.')
                WriteFile = True
            else:
                if not quiet:
                    print('>>> Will not overwrite.')
                WriteFile = False

        # if file doesn't exist, allow to write
        else:
            WriteFile = True

        return WriteFile



    
    # check inputs and check file exists
    #----------------------------------
    if not isinstance(hemisphere, str):
        raise TypeError(f'hemisphere should be string, not {type(hemisphere)}')
    if hemisphere not in ['N', 'S']:
        raise ValueError(f"hemisphere should be 'N' or 'S', not '{hemisphere}'")

    # open specified file
    if not quiet: 
        print(f'>>> opening {file}')
    
    # check whether file exists
    file_exists = True
    if not os.path.exists(file):
        file_exists = False
    assert file_exists, f'!!!  {file} does not exist'
            
        
    # open data
    #----------------------------------
    f = h5py.File(file,'r')

    # navigate to data of chosen hemisphere
    F = f['HDFEOS']['GRIDS'][f'{hemisphere}pPolarGrid06km']

    # grab geo data
    #---------------
    lat = F['lat'][:]
    lon = F['lon'][:]
    lon[lon<0]+=360

    # find units
    lat_units = F['lat'].attrs['units'].decode()
    lon_units = F['lon'].attrs['units'].decode()
    

    # grab specified Tb data
    #-----------------------
    data_vars = ['89H_DAY', '89V_DAY', '89H_ASC', '89V_ASC', '89H_DSC','89V_DSC']
    vars_dict = {}

    for VAR in data_vars:

        vars_dict[VAR] = {}

        # grab data and packing convention info
        data = F['Data Fields'][f'SI_06km_NH_{VAR}'][:].astype(float)
        _FillValue   = F['Data Fields'][f'SI_06km_NH_{VAR}'].attrs['_FillValue'][0]
        scale_factor = F['Data Fields'][f'SI_06km_NH_{VAR}'].attrs['scale_factor'][0]
        add_offset   = F['Data Fields'][f'SI_06km_NH_{VAR}'].attrs['add_offset'][0]

        # unpack data
        data[data == _FillValue] = np.nan
        scaled_data = ( scale_factor * data ) + add_offset

        # grab data units and variable name
        data_units = F['Data Fields'][f'SI_06km_NH_{VAR}'].attrs['units'].decode()
        data_name = F['Data Fields'][f'SI_06km_NH_{VAR}'].attrs['long_name'].decode()

        # add to dict
        vars_dict[VAR]['data'] = scaled_data
        vars_dict[VAR]['units'] = data_units
        vars_dict[VAR]['long_name'] = data_name

        
    # crop all data to specified range
    #---------------------------------
    if not quiet: 
        print(f'>>> cropping lat: {lat_range[0]}-{lat_range[1]} lon: {lon_range[0]}-{lon_range[1]}')
    lon, lat, VARS = crop_data(lon=lon, lat=lat, 
                               VARS=[vars_dict['89H_DAY']['data'], vars_dict['89V_DAY']['data'],
                                     vars_dict['89H_ASC']['data'], vars_dict['89V_ASC']['data'],
                                     vars_dict['89H_DSC']['data'], vars_dict['89V_DSC']['data']], 
                               lon_range=lon_range, lat_range=lat_range)
        
    # grab file name, path, and date
    #---------------------------------
    # separate file path and name
    filename = file.split('/')[-1]
    filepath = file[:file.find(filename)]

    # grab file date string from file name
    date_string = filename.split('.')[0].split('_')[-1]
    date = datetime.datetime.strptime(date_string, '%Y%m%d')
    
    
    # create netCDF data set
    #-----------------------
    ds = xr.Dataset(
    data_vars=dict(
        DAY_89H = (["x", "y"], VARS[0], 
                   {"description":vars_dict['89H_DAY']['long_name'], "units":vars_dict['89H_DAY']['units']}),
        DAY_89V = (["x", "y"], VARS[1], 
                   {"description":vars_dict['89V_DAY']['long_name'], "units":vars_dict['89V_DAY']['units']}),
        ASC_89H = (["x", "y"], VARS[2], 
                   {"description":vars_dict['89H_ASC']['long_name'], "units":vars_dict['89H_ASC']['units']}),
        ASC_89V = (["x", "y"], VARS[3], 
                   {"description":vars_dict['89V_ASC']['long_name'], "units":vars_dict['89V_ASC']['units']}),
        DSC_89H = (["x", "y"], VARS[4], 
                   {"description":vars_dict['89H_DSC']['long_name'], "units":vars_dict['89H_DSC']['units']}),
        DSC_89V = (["x", "y"], VARS[5], 
                   {"description":vars_dict['89V_DSC']['long_name'], "units":vars_dict['89V_DSC']['units']}),
    ),
    coords=dict(
        longitude = (["x", "y"], lon, {"units": lon_units}),
        latitude  = (["x", "y"], lat, {"units": lat_units}),
    ),
    attrs=dict(
        description = (f"6.25 km 89 GHz brightness temperatures, cropped ([{lat_range[0]},{lat_range[1]}], [{lon_range[0]},{lon_range[1]}) and converted from large (122 MB) locally-stored .he5 NSIDC data (NSIDC-AU_SI6, doi: 10.5067/NX1R09ORNOZN)"), 
        time = (f"{date}")),
)
        
        
    # create name of new file and locate path to save it to
    #------------------------------------------------------
    # create full name of file to save
    savename = new_filename.format(date_string)
    out_file = filepath+save_folder+'/'+savename

    # create savefolder if it does not already exist
    if not os.path.exists(filepath+save_folder):
        os.mkdir(filepath+save_folder)
        print(f'>>> created new folder to save .nc files: {filepath+save_folder}')


    # save file
    #------------------------------------
    if not quiet:
        print(f'>>> saving file: {savename}')
        print(f'>>> saving to: {filepath+save_folder+"/"}')

    # check whether the file already exists, and if so, whether can be overwritten
    WriteFile = check_writefile(out_file, allow_overwrites = allow_overwrites, quiet = quiet)

    # write file if allowed
    if WriteFile:
        if not quiet:
            print(f'\nSaving {out_file}')
        ds.to_netcdf(path=out_file, mode='w', format="NETCDF4")

    # move he5 file to save_folder if specified
    #------------------------------------------
    if move_he5_file:

        # create  name of .he5 file once moved to destination
        moved_he5_file = filepath+save_folder+'/'+filename

        # check whether the file already exists, and if so, whether can be overwritten
        WriteFile = check_writefile(moved_he5_file, allow_overwrites = allow_overwrites, quiet = quiet)

        # write new file
        #---------------
        if WriteFile:
            # move file
            os.rename(file, moved_he5_file)
            if not quiet:
                print(f'>>> moved "{filename}"\n    "{filepath}" --> "{filepath+save_folder+"/"}"')

    # return nc file name if specified
    if return_new_filename:
        return out_file
    
    







#////////////////////
#  read_AMSR_he5 ///
#//////////////////
#---------------------------------------------------------------------
# Grab 89 GHz brightness temperatures from locally-stored .he5 data
#---------------------------------------------------------------------
# DEPENDENCIES:
import os
import h5py
import datetime
from metpy.units import units
import numpy as np
#---------------------------------------------------------------------
def read_AMSR_he5(date = datetime.datetime(year = 2000, month = 1, day = 1),
                  filepath = '/Users/mackenziejewell/Downloads/', 
                  filename_convention = 'AMSR_U2_L3_SeaIce6km_B04_{}{}{}.he5', 
                  hemisphere = 'N', add_units = False,
                  return_vars = ['lat', 'lon', '89H_DAY', '89V_DAY', '89H_ASC', '89V_ASC', '89H_DSC', '89V_DSC'],    
                  quiet = True):
    """Grab 6.25 km 89 GHz brightness temperatures from locally-stored .he5 NSIDC data (NSIDC-AU_SI6, doi: 10.5067/NX1R09ORNOZN)
INPUT: 
- date: date of data to open
- filepath: path to local .he5 data files
- filename_convention: naming convention of stored data files, 
    with {}{}{} indicating in order year / month / day as YYYY / MM / DD
- hemisphere: either 'N' for north hemisphere data or 'S' for south hemisphere data
- add_units: whether or not to add metpy units to data (default: False)
    * NOTE! this is much slower to display. Can display variable quick with variable.magnitude
- return_vars: variables/attributes to return in specified order. (list)
    Can include any or all OUTPUT variables in any order: 
    default: ['lat', 'lon', '89H_DAY', '89V_DAY', '89H_ASC', '89V_ASC', '89H_DSC', '89V_DSC']

- quiet: bool, whether or not to hide print statements (default: True)
    

OUTPUT:
List of any or all of variables specified in return_vars:
(or list of empty lists if data file does not exist)
- lat: (MxN) grid of latitude data
- lon: (MxN) grid of longitude data (0-360)
- 89H_DAY: (MxN) grid of H pol 89 GHz daily average Tb
- 89V_DAY: (MxN) grid of V pol 89 GHz daily average Tb
- 89H_ASC: (MxN) grid of H pol 89 GHz daily average ascending Tb
- 89V_ASC: (MxN) grid of V pol 89 GHz daily average ascending Tb
- 89H_DSC: (MxN) grid of H pol 89 GHz daily average descending Tb
- 89V_DSC: (MxN) grid of V pol 89 GHz daily average descending Tb


DEPENDENCIES:
import os
import h5py
import datetime
from metpy.units import units
import numpy as np


Latest recorded update:
08-02-2023
    """
    
    
    # check inputs and check file exists
    #----------------------------------
    if not isinstance(hemisphere, str):
        raise TypeError(f'hemisphere should be string, not {type(hemisphere)}')
    if hemisphere not in ['N', 'S']:
        raise ValueError(f"hemisphere should be 'N' or 'S', not '{hemisphere}'")
    
    if not quiet: 
        print(f'>>> open file for date: {date}')
    
    # open file corresponding to date
    filename = filename_convention.format(date.year, str(date.month).zfill(2), str(date.day).zfill(2))
    
    if not quiet: 
        print(f'>>> opening {filename}')
    
    # create full file path + name
    file = filepath + filename
    
    # check whether file exists
    file_exists = True
    if not os.path.exists(file):
        file_exists = False
        if not quiet:
            print(f'!!!  {file} does not exist')

            
    # open data
    #----------------------------------
    
    # store all output in dictionary
    vars_dict = {}
    
    
    # if file exists, open relevant data
    if file_exists:
        f = h5py.File(file,'r')
        
        # navigate to data of chosen hemisphere
        F = f['HDFEOS']['GRIDS'][f'{hemisphere}pPolarGrid06km']
        
        
        # grab geo data if specified
        
        # latitude
        #---------
        if 'lat' in return_vars:
            # grab from data
            lat = F['lat'][:]
            # add units if specified
            if add_units:
                lat = lat * units(F['lat'].attrs['units'].decode())
            # add to dict
            vars_dict['lat'] = lat
            
        # longitude
        #----------
        if 'lon' in return_vars:
            # grab from data
            lon = F['lon'][:]
            lon[lon<0]+=360
            # add units if specified
            if add_units:
                lon = lon * units(F['lon'].attrs['units'].decode())
            # add to dict
            vars_dict['lon'] = lon
                
                
        # grab specified Tb data
        #-----------------------
        for VAR in return_vars:
            if VAR not in ['lat', 'lon']:
                # grab data and packing convention info
                data = F['Data Fields'][f'SI_06km_NH_{VAR}'][:].astype(float)
                _FillValue   = F['Data Fields'][f'SI_06km_NH_{VAR}'].attrs['_FillValue'][0]
                scale_factor = F['Data Fields'][f'SI_06km_NH_{VAR}'].attrs['scale_factor'][0]
                add_offset   = F['Data Fields'][f'SI_06km_NH_{VAR}'].attrs['add_offset'][0]

                # unpack data
                data[data == _FillValue] = np.nan
                scaled_data = ( scale_factor * data ) + add_offset

                # add units if specified
                if add_units:
                    scaled_data = scaled_data * units(F['Data Fields'][f'SI_06km_NH_{VAR}'].attrs['units'].decode())

                # add to dict
                vars_dict[VAR] = scaled_data

        
    # if file does not exist, return all vars as empty lists
    else:
        vars_dict['lon'] = []
        vars_dict['lat'] = []
        vars_dict['89H_DAY'] = []
        vars_dict['89V_DAY'] = []
        vars_dict['89H_ASC'] = []
        vars_dict['89V_ASC'] = []
        vars_dict['89H_DSC'] = []
        vars_dict['89V_DSC'] = []

    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    if len(return_data) == 1:
        return_data = return_data[0]
    return return_data



#////////////////
#  crop_data ///
#//////////////
#---------------------------------------------------------------------
# Crop AMSR data and coords to within given lat/lon range
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def crop_data(lon=[], lat=[], VARS=[], lon_range=[180, 260], lat_range=[65,90]):
    
    """Crop lats, lons, variables on non-equirectangular grid to within given lat/lon range.
    
INPUT:
- VARS: (Lx1) list of variables, each with shape (MxN). Even if just one variable, enclose in brackets as list.
- lon: (MxN) longitude grid (0 to 360) associated with VARS list
- lat: (MxN) latitude grid associated with VARS lists
- lon_range:longitude range for cropping (defined 0 to 360)
    (default: [180, 260])
- lat_range:latitude range for cropping
    (default: [65,90])

OUTPUT:
- lon, lat, VARS  --> cropped input grids
    
DEPENDENCIES:
import numpy as np
import numpy.ma as ma

Latest recorded update:
08-03-2023
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
        lat = np.delete(lat, rows_to_delete, 0)
        lon = np.delete(lon, rows_to_delete, 0)
        for vv in range(len(VARS)):
            VARS[vv] = np.delete(VARS[vv], rows_to_delete, 0)

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
        lat = np.delete(lat, columns_to_delete, 1)
        lon = np.delete(lon, columns_to_delete, 1)
        for vv in range(len(VARS)):
            VARS[vv] = np.delete(VARS[vv], columns_to_delete, 1)

        # determine change in array sizes to see if it is still cropping
        array_size_change = (num_rows_i-np.shape(lat)[0])+(num_cols_i-np.shape(lat)[1])

    return lon, lat, VARS

        