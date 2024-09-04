from pyhdf.SD import SD, SDC
import re
import numpy as np
import glob
import cartopy.crs as ccrs

def open_MODIS_IST_file(file):
    
    # Read dataset.
    hdf = SD(file, SDC.READ)
    IST_raw = hdf.select(f'Ice_Surface_Temperature')[:].astype(float)
    IST_QA = hdf.select(f'Ice_Surface_Temperature_Spatial_QA')[:].astype(float)
    
    scale_factor = 0.01
    add_offset = 0
    IST = scale_factor * (IST_raw - add_offset)
    valid_range = (210, 313.2)
    IST[IST<valid_range[0]] = np.nan
    IST[IST>valid_range[1]] = np.nan

    land_mask = (IST_QA == 253)
    ocean_mask = (IST_QA == 254)
    bad_quality = (IST_QA == 1)
    
    # Construct the grid.  The needed information is in a global attribute
    # called 'StructMetadata.0'.  Use regular expressions to tease out the
    # extents of the grid.
    fattrs = hdf.attributes(full=1)
    ga = fattrs["StructMetadata.0"]
    gridmeta = ga[0]
    ul_regex = re.compile(r'''UpperLeftPointMtrs=\(
                              (?P<upper_left_x>[+-]?\d+\.\d+)
                              ,
                              (?P<upper_left_y>[+-]?\d+\.\d+)
                              \)''', re.VERBOSE)

    match = ul_regex.search(gridmeta)
    x0 = float(match.group('upper_left_x'))
    y0 = float(match.group('upper_left_y'))

    lr_regex = re.compile(r'''LowerRightMtrs=\(
                              (?P<lower_right_x>[+-]?\d+\.\d+)
                              ,
                              (?P<lower_right_y>[+-]?\d+\.\d+)
                              \)''', re.VERBOSE)
    match = lr_regex.search(gridmeta)
    x1 = float(match.group('lower_right_x'))
    y1 = float(match.group('lower_right_y'))

    nx, ny = IST_raw.shape
    x = np.linspace(x0, x1, nx, endpoint=False)
    y = np.linspace(y0, y1, ny, endpoint=False)
    xv, yv = np.meshgrid(x, y)
    
    return x, y, IST

def find_MODIS_IST_files(date, file_path = ''):
    
    # find all hdf files
    hdf_files = glob.glob(file_path+'*.hdf')

    files = []
    tiles = []

    # loop throught files, select only those matching 
    # 'MYD29P1N' group and desired date
    for file in hdf_files:

        if 'MYD29P1N' in file:
            file_date = file[file.find('.A')+2:].split('.')[0]

            if file_date == date.strftime('%Y%j'):
                file_tile = file[file.find('.A')+2:].split('.')[1]
                files.append(file)
                tiles.append(file_tile)
            
    return files, tiles
    
    

def grab_MODIS_IST(date, file_path = ''):
    
    files, tiles = find_MODIS_IST_files(date, file_path = file_path)

    tile_order = ['h07v07', 'h08v07']

    # open first file
    X1, Y1, IST1 = open_MODIS_IST_file(files[tiles.index('h07v07')])
    X_grid1, Y_grid1 = np.meshgrid(X1, Y1)

    # open second file
    X2, Y2, IST2 = open_MODIS_IST_file(files[tiles.index('h08v07')])
    X_grid2, Y_grid2 = np.meshgrid(X2, Y2)

    # stack grids
    Xgrid = np.concatenate((X_grid1,X_grid2),axis=1)
    Ygrid = np.concatenate((Y_grid1,Y_grid2),axis=1)
    IST = np.concatenate((IST1,IST2),axis=1)

    # create ice projection from info (EASE-grid)
    # based off documentation https://nsidc.org/sites/default/files/myd29p1n-v061-userguide_1.pdf
    ice_projection = ccrs.LambertAzimuthalEqualArea(central_longitude=0, central_latitude=90,
                                                    globe=ccrs.Globe(semimajor_axis = 6371228, semiminor_axis = 6371228))

    data = {}
    data['xx'] = Xgrid
    data['yy'] = Ygrid
    data['ist'] = IST
    data['proj'] = ice_projection
    
    return data
    