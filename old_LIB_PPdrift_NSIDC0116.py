# OLD versions of routines used to open PP sea ice drift data

#////////////////////
#  get_PPD_date  ///
#//////////////////
#---------------------------------------------------------------------
# Find date index of given day in given NSIDC Polar Pathfinder file. 
#---------------------------------------------------------------------
#////////////////////
#  get_PPD_data  ///
#//////////////////
#---------------------------------------------------------------------
# Grab the u/v (polar EASE grid components) data and associated lat/lon grids 
# at the time associated with PPD_index from NSIDC Polar Pathfinder file
#---------------------------------------------------------------------
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

#////////////////////
#  get_PPD_date  ///
#//////////////////
#---------------------------------------------------------------------
# Find date index of given day in given NSIDC Polar Pathfinder file. 
#---------------------------------------------------------------------
# DEPENDENCIES:
#-------------
import netCDF4
#---------------------------------------------------------------------
def get_PPD_date(dt_obj, PPD_file):
    
    """Find date index of given day in given NSIDC Polar Pathfinder file. Given datetime object ,return nearest date as cftime object (and associated index) for desired date in file. Throw error and return None, None if no day match is found.
    
INPUT:
- dt_obj: datetime object of desired date to be found in PP file
- PPD_file: PolarPathfinder daily drift netCDF file (with path)

OUTPUT:
- PPD_date: nearest date to dt_obj in PP data (as cftime.DatetimeJulian object)
- PPD_index: index of PP_date in PP data file

DEPENDENCIES:
import netCDF4

Latest recorded update:
05-03-2022
    """
    
    # grab times from ice_file
    #-------------------------
    ice_file = netCDF4.Dataset(PPD_file, 'r')
    time = ice_file.variables['time']
    dates = netCDF4.num2date(time[:], time.units, time.calendar)

    # first set PP_index as None in case day match cannot be found in file
    #---------------------------------------------------------------------
    PPD_index = None
    # check for a match in the year, month, and day of Polar Pathfinder file
    #-----------------------------------------------------------------------
    for d in range(len(dates)):
        if dt_obj.day == dates[d].day and dt_obj.month == dates[d].month and dt_obj.year == dates[d].year:
            PPD_index = d
            break
    
    # if no day match found, throw error
    # and return None
    #-----------------------------------
    if PPD_index == None:
        print('ERROR: No day match to {} found in {}'.format(dt_obj, PPD_file))
        return None, None
    
    # else grab and return date
    #--------------------------
    else:
        PPD_date = dates[PPD_index]
        return PPD_date, PPD_index
    
    # close file
    #-----------
    ice_file.close()
    
    
    
    
    
    
    
#////////////////////
#  get_PPD_data  ///
#//////////////////
#---------------------------------------------------------------------
# Grab the u/v (polar EASE grid components) data and associated lat/lon grids 
# at the time associated with PPD_index from NSIDC Polar Pathfinder file
#---------------------------------------------------------------------
def get_PPD_data(PPD_file, PPD_index):
    
    """Grab the u/v (polar EASE grid components) data and associated lat/lon grids 
    at the time associated with PPD_index from NSIDC Polar Pathfinder file. ref: https://nsidc.org/support/how/how-convert-horizontal-and-vertical-components-east-and-north to understand how the u/v components relate to eastward/northward.
    
INPUT:
- PPD_file: PolarPathfinder daily drift netCDF file (with path)
- PP_index: index of desried date in PP data file

OUTPUT:
- u, v: velocity components on EASE grid (cm/s)
    (u: toward the right on the grid)
    (v: upward (toward the top) on the grid)
- lons, lats: spatial grids associated with u, v
    lons are from (0 to 360)
    

DEPENDENCIES:
import netCDF4
import numpy as np

Latest recorded update:
05-03-2022
    """
    
    # open PP file
    #-------------
    ice_file = netCDF4.Dataset(PPD_file, 'r')
    
    # grab data at PP_index time
    #---------------------------
    u = ice_file.variables['u'][PPD_index,:,:]
    v = ice_file.variables['v'][PPD_index,:,:]
    lats = ice_file.variables['latitude'][:].data
    lons = ice_file.variables['longitude'][:].data
    
#     # check that no ice velocities are outside 
#     # the valid_min or valic_max range
#     #-----------------------------------------
#     u_min,u_max = ice_file.variables['u'].valid_min, ice_file.variables['u'].valid_max
#     v_min,v_max = ice_file.variables['v'].valid_min, ice_file.variables['v'].valid_max
#     for ii in range(np.shape(u)[0]):
#         for jj in range(np.shape(u)[1]):
#             if u[ii][jj] < u_min or u[ii][jj] > u_max:
#                 u[ii][jj] = np.nan
#             if v[ii][jj] < v_min or v[ii][jj] > v_max:
#                 v[ii][jj] = np.nan

    # make all longitudes range (0,360)
    #----------------------------------
    for i in range(lons.shape[0]):
        for j in range(lons.shape[0]):
            if lons[i][j] < 0:
                lons[i][j] += 360                 
    
    # close file
    #-----------
    ice_file.close()
    
    return lons, lats, u, v



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
def convert_PPD_vectors(lons, u, v, fill_to_nan = False):
    
    """Convert u/v (polar EASE grid) velocity components from NSIDC Polar Pathfinder data to 
    eastward and northward vector components. ref: https://nsidc.org/support/how/how-convert-horizontal-and-vertical-components-east-and-north.
    
INPUT:
- u, v: velocity components on EASE grid (cm/s)
    (u: toward the right on the grid)
    (v: upward (toward the top) on the grid)
- lons: longitude grid (0 to 360) associated with u, v
- fill_to_nan: specify whether or not to replace fill value in masked arrays 
    with nans in calculation so fill number is not treated as actual data value. 
    (default: False if no changes applied to data)
    (else set to True and routine will determine masked array fill value and replace with nans).

OUTPUT:
- E: grid of east component of velocity (cm/s)
- N: grid of north component of velocity (cm/s)
    

DEPENDENCIES:
import numpy as np

Latest recorded update:
05-03-2022
    """
                
    if fill_to_nan == True:
        u_fill = u.fill_value
        u = u.data
        u[u==u_fill] = np.nan
        
        v_fill = v.fill_value
        v = v.data
        v[v==fill_to_nan] = np.nan
        
    # convert EASE grid vector components to northward, eastward
    #------------------------------------------------------------------
    E = u * np.cos(lons/180*np.pi)  +  v * np.sin(lons/180*np.pi)
    N = -u * np.sin(lons/180*np.pi)  +  v * np.cos(lons/180*np.pi)
    
                
    return E, N

#////////////////////
#  crop_PPD_data ///
#//////////////////
#---------------------------------------------------------------------
# Crop NSIDC Polar Pathfinder lats, lons, u, v, to within given lat/lon range
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def crop_PPD_data(lons, lats, u, v, lon_range=[195, 235], lat_range=[65,78]):
    
    """Crop NSIDC Polar Pathfinder lats, lons, u, v, to within given lat/lon range.
    
INPUT:
- u: grid of eastward vector components of ice drift
- v: grid of northward vector components of ice drift
- lons: longitude grid (0 to 360) associated with u, v
- lats: latitude grid associated with u, v
- lon_range:longitude range for cropping (defined 0 to 360)
    (default: [195, 235])
- lat_range:latitude range for cropping
    (default: [65,78])

OUTPUT:
- lons, lats, u, v  --> cropped input grids
    
DEPENDENCIES:
import numpy as np
import numpy.ma as ma

Latest recorded update:
05-03-2022
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
        num_rows_i = np.shape(lats)[0]
        num_cols_i = np.shape(lats)[1]
        # run through each row in lat/lons and check 
        # if any coords are within desired range, keep row
        #--------------------------------------------------
        rows_to_delete = []
        for ii in range(np.shape(lats)[0]):
            # check whether all lats/lons in row are within ranges
            row_check = 0
            for jj in range(np.shape(lats)[1]):
                # add 1 to row_check and break current loop if any coords are within ranges
                if lats[ii][jj] < latmax and lats[ii][jj] > latmin and lons[ii][jj] < lonmax and lons[ii][jj] > lonmin:
                    row_check+=1
                    break
            # if no coords values within desired range were found in row,
            # add row index to rows_to_delete
            if row_check == 0:
                rows_to_delete.append(ii)

        # crop lats, lons, u, v to new range
        #------------------------------------
        lats = np.delete(lats, rows_to_delete, 0)
        lons = np.delete(lons, rows_to_delete, 0)
        u = np.delete(u, rows_to_delete, 0)
        v = np.delete(v, rows_to_delete, 0)

        # run through each column in lat/lons and check 
        # if any coords are within desired range, keep column
        #----------------------------------------------------
        columns_to_delete = []
        for jj in range(np.shape(lats)[1]):
            # check whether all lats/lons in column are within ranges
            column_check = 0
            for ii in range(np.shape(lats)[0]):
                # add 1 to column_check and break current loop if any
                # coords are within ranges
                if lats[ii][jj] < latmax and lats[ii][jj] > latmin and lons[ii][jj] < lonmax and lons[ii][jj] > lonmin:
                    column_check+=1
                    break
            # if no coords values within desired range were found in column, add column index to columns_to_delete
            if column_check == 0:
                columns_to_delete.append(jj)

        # crop lats, lons, u, v to new range
        #------------------------------------
        lats = np.delete(lats, columns_to_delete, 1)
        lons = np.delete(lons, columns_to_delete, 1)
        u = np.delete(u, columns_to_delete, 1)
        v = np.delete(v, columns_to_delete, 1)

        # determine change in array sizes to see if it is still cropping
        array_size_change = (num_rows_i-np.shape(lats)[0])+(num_cols_i-np.shape(lats)[1])

    # np.delete messes with mask, so remask
    u = ma.masked_where(u == -9999, u)
    v = ma.masked_where(v == -9999, v)
    
    return lons, lats, u, v



#/////////////////////
#  grab_ice_Drift ///
#///////////////////
#---------------------------------------------------------------------
# Import NSIDC Polar Pathfinder lats, lons, u, v cropped to within given lat/lon range.
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def grab_ice_Drift(dt_obj, PPD_drift_path = '/Users/mackenziejewell/Data/PP_icedrift/', 
                   PPD_filename = 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc', 
                  lat_range = [65,78], lon_range = [195, 235]):
    
    """Import NSIDC Polar Pathfinder lats, lons, u, v cropped to within given lat/lon range.
    
INPUT:
- dt_obj: desired date
- PPD_drift_path: 
- PPD_filename: naming convention for PPD files (default: 'icemotion_daily_nh_25km_{}0101_{}1231_v4.1.nc' where {} will be replace with year of dt_obj)
- lat_range:latitude range for cropping
    (default: [65,78])
- lon_range:longitude range for cropping (defined 0 to 360)
    (default: [195, 235])

OUTPUT:
- lons: M x N longitude grid (0 to 360) associated with u, v
- lats: M x N latitude grid associated with u, v
- u: M x N grid of eastward vector components of ice drift
- v: M x N grid of northward vector components of ice drift
    
DEPENDENCIES:
import numpy as np
import numpy.ma as ma

Latest recorded update:
06-03-2022
    """
    PP_mainpath = PPD_drift_path+PPD_filename
    PP_file = PP_mainpath.format(dt_obj.year,dt_obj.year)
    PP_date, PP_index = get_PPD_date(dt_obj, PP_file)
    lons, lats, u, v = get_PPD_data(PP_file, PP_index)
    lons, lats, u, v = crop_PPD_data(lons, lats, u, v, lon_range=lon_range, lat_range=lat_range)
    u, v = convert_PPD_vectors(lons, u, v, fill_to_nan = True)
    
    return lons, lats, u, v
                                    