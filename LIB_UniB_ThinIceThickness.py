# DEPENDENCIES:
import os
from datetime import datetime
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
import xarray as xr

def grab_projinfo_tit(ds, suppress_prints = True):

    """Grab projection info from smos-smap thin ice thickness data (doi: 10.5194/tc-13-675-2019)
    https://seaice.uni-bremen.de/thin-ice-thickness/

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
08-22-2024
    """
    
    
    # grab projection info
    #---------------------
    CRS = ds.crs.attrs

    # grab parameters from crs spatial attributes
    central_meridian = int(CRS['straight_vertical_longitude_from_pole'])
    semimajor = CRS['semi_major_axis']
    semiminor = CRS['semi_minor_axis']
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
                                           globe=ccrs.Globe(semimajor_axis = semimajor,
                                                            semiminor_axis = semiminor, inverse_flattening=inv_flat),
                                           true_scale_latitude=standard_parallel)

    return projection



def grab_ThinIceThickness(date = datetime(2018, 1, 1), 
                          main_path = '/Volumes/Seagate_Jewell/KenzieStuff/UniB-SMOS-SMAP/', 
                          file_convention = '{}_north_mix_sit_v300.nc'):
    
    # construct filename
    #-------------------
    # annual folder
    year_folder = main_path + f'{date.year}/'
    
    # convert date to string for filename
    date_string = date.strftime('%Y%m%d')
    filename = file_convention.format(date_string)
    
    # open with xarray 
    #-------------------
    ds = xr.open_dataset(year_folder + filename)
    ds.close()
    
    # read data variables into dictionary
    #-------------------
    data = {}
    data['x'] = ds.x.values
    data['y'] = ds.y.values
    data['xx'], data['yy'] = np.meshgrid(data['x'], data['y'])

    for var in ['smos_thickness', 'smos_thickness_unc', 'smap_thickness', 'smap_thickness_unc']:
        data[var] = ds[var].values

    # combined thickness
    data['tit'] =  ds['combined_thickness'].values
    data['tit_unc'] =  ds['combined_thickness_unc'].values

    data['proj'] = grab_projinfo_tit(ds, suppress_prints=True)
    
    data['ds'] = ds
    
    return data