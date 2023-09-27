#/////////////////////////////
#  convert_AMSR_hdf_to_nc ///
#///////////////////////////
#---------------------------------------------------------------------
# Crop and convert locally-stored .hdf files to smaller .nc files
#---------------------------------------------------------------------
#////////////////
#  crop_data ///
#//////////////
#---------------------------------------------------------------------
# Crop AMSR data and coords to within given lat/lon range
#---------------------------------------------------------------------


#/////////////////////////////
#  convert_AMSR_hdf_to_nc ///
#///////////////////////////
#---------------------------------------------------------------------
# Crop and convert locally-stored .hdf files to smaller .nc files
#---------------------------------------------------------------------
# DEPENDENCIES:
import os
from pyhdf.SD import SD, SDC
import numpy as np
import xarray as xr
import datetime
#---------------------------------------------------------------------
def convert_AMSR_hdf_to_nc(file = '/Volumes/Seagate_Jewell/KenzieStuff/AMSR89Tb_NSIDCAE_SI6/AMSR_E_L3_SeaIce6km_V15_20030101.hdf',
                           AMSR_latlon_file = '/Volumes/Seagate_Jewell/KenzieStuff/AMSR89Tb_NSIDCAE_SI6/NSIDC0771_LatLon_PS_N6.25km_v1.0.nc',
                           hemisphere = 'N', lon_range=[180, 260], lat_range=[65,90], 
                           new_filename = 'AMSR_89Ghz_{}.nc', save_folder = 'converted',
                           move_hdf_file = False, allow_overwrites = False,
                           return_new_filename = False,
                           quiet = True):
    """Grab 6.25 km 89 GHz brightness temperatures from locally-stored AMSR-E NSIDC data (NSIDC-AE_SI6, doi: 10.5067/AMSR-E/AE_SI6.003) and crop then save as much smaller .nc file.
    
INPUT: 
- file: name of stored .hdf data file, including path to file
    (default: '/Volumes/Seagate_Jewell/KenzieStuff/AMSR89Tb_NSIDCAE_SI6/AMSR_E_L3_SeaIce6km_V15_20030101.hdf')
- AMSR_latlon_file: name of stored .hdf data file containing lat/lon coords of NSIDC NPS projection, including path to file
    (default: '/Volumes/Seagate_Jewell/KenzieStuff/AMSR89Tb_NSIDCAE_SI6/NSIDC0771_LatLon_PS_N6.25km_v1.0.nc')
- hemisphere: string specifying hemispheric data to grab.
    either 'N' for north hemisphere data or 'S' for south hemisphere data
- lon_range: [lonmin, lonmax] to crop, values in range (0,360). default: [180, 260]
- lat_range: [latmin, latmax] to crop. default: [65,90]

- new_filename: name with which to save new file (with {} inserted where date info will be placed)
        (default: 'AMSR_89Ghz_{}.nc')
- save_folder: name of folder within file's path to save file to. Will then copy used .he5 file to same folder.
                Do not add any slashes to folder name. (default: 'converted').
                If folder does not yet exist, will create folder.
- move_hdf_file: bool, whether or not to move .hdf file to save_folder alongside new .nc file
    (default: False)              
- allow_overwrites: bool, whether or not to allow file overwrites (default: False)
- return_new_filename: bool, whether or not to return name of new .nc file. 
    if True, return filename. 
    if False (default), don't return anything.
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
from pyhdf.SD import SD, SDC
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
    try:
        f = SD(file,SDC.READ) 
    # raise an error if file can't be opened 
    except Exception as e:
        print(e, ", error opening ", file, sep='')

    # grab geo data
    #---------------
    # use polar stereographic coordinate file to grab data coords
    # https://nsidc.org/data/nsidc-0771/versions/1
    # NSIDC0771_LatLon_PS_N6.25km_v1.0.nc file 
    ds = xr.open_dataset(AMSR_latlon_file)

    # grab lats/lons of projected data
    lat = ds.latitude.values
    lon = ds.longitude.values
    lon[lon<0]+=360

    # grab x, y if desired
    xx = ds.x.values
    yy = ds.y.values
    CRS = ds.crs.attrs # projection dict

    # define units
    lat_units = 'degrees_north'
    lon_units = 'degrees_east'
    
    
    # grab specified Tb data
    #-----------------------
    data_vars = ['89H_DAY', '89V_DAY', '89H_ASC', '89V_ASC', '89H_DSC','89V_DSC']
    vars_dict = {}

    for VAR in data_vars:
        
        vars_dict[VAR] = {}

        # grab data and packing convention info
        data = f.select(f'SI_06km_{hemisphere}H_{VAR}')[:].astype(float)
        
        # according to documentation:
        # https://nsidc.org/sites/default/files/ae_si6-v003-userguide.pdf
        _FillValue   = 0
        scale_factor = 0.1
        add_offset   = 0

        # unpack data
        data[data == _FillValue] = np.nan
        scaled_data = ( scale_factor * data ) + add_offset

        # add to dict
        #--------------------------------
        vars_dict[VAR]['data'] = scaled_data
        vars_dict[VAR]['units'] = 'degree_kelvin'
        
        # ascending, descending, or daily
        if "ASC" in VAR:
            des = 'ascending '
        elif "DSC" in VAR:
            des = 'descending '
        else:
            des = ''
            
        # vertical or horizontal polarize
        if "V" in VAR:
            pol = 'vertical'
        else:
            pol = 'horizontal'
            
        vars_dict[VAR]['long_name'] = f'89.0 GHz {pol} daily average {des}Tbs'
        
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
        description = (f"6.25 km 89 GHz brightness temperatures, cropped ([{lat_range[0]},{lat_range[1]}], [{lon_range[0]},{lon_range[1]}) and converted from large (47 MB) locally-stored .hdf NSIDC data (NSIDC-AE_SI6, doi: 10.5067/AMSR-E/AE_SI6.003)"), 
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

    # move hdf file to save_folder if specified
    #------------------------------------------
    if move_hdf_file:

        # create  name of .hdf file once moved to destination
        moved_hdf_file = filepath+save_folder+'/'+filename

        # check whether the file already exists, and if so, whether can be overwritten
        WriteFile = check_writefile(moved_hdf_file, allow_overwrites = allow_overwrites, quiet = quiet)

        # write new file
        #---------------
        if WriteFile:
            # move file
            os.rename(file, moved_hdf_file)
            if not quiet:
                print(f'>>> moved "{filename}"\n    "{filepath}" --> "{filepath+save_folder+"/"}"')

    # return nc file name if specified
    if return_new_filename:
        return out_file



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

        