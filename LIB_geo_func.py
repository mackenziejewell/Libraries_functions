#////////////////////
#  make_polygon  ///
#//////////////////
#---------------------------------------------------------------------
# create shapely polygon from coordinates of desired polygon perimeter
#---------------------------------------------------------------------
#///////////////////////////////////
#  make_smooth_geographic_box   ///
#/////////////////////////////////
#---------------------------------------------------------------------
# create smooth geographic box that smoothly follows parallels and 
# meridians between input coordinates
#---------------------------------------------------------------------
#///////////////////////
#  make_linestring  ///
#/////////////////////
#---------------------------------------------------------------------
# create shapely linsetring from array of coordinates
#---------------------------------------------------------------------
#///////////////////////
#  regrid_geo_data  ///
#/////////////////////
#---------------------------------------------------------------------
# regrid data given initial data grid and desired geo grid
#---------------------------------------------------------------------
#//////////////////////////////
#  within_polygon_indices  ///
#////////////////////////////
#---------------------------------------------------------------------
# Find indices of geo grid that fall within polygon.
#---------------------------------------------------------------------
#/////////////////////////////////
#  polygon_mean_from_indices  ///
#///////////////////////////////
#---------------------------------------------------------------------
# Find mean value of data within polygon from indices of data_grid
#---------------------------------------------------------------------




#////////////////////
#  make_polygon  ///
#//////////////////
#---------------------------------------------------------------------
# create shapely polygon from coordinates of desired polygon perimeter
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np, numpy.ma as ma
from pyproj import Geod
from shapely import wkt
#---------------------------------------------------------------------

def make_polygon(perim_coords, ellipsoid = "WGS84", quiet = True):
    
    """Create shapely polygon from coordinates of desired polygon perimeter.
    
INPUT:
- perim_coords: Nx2 array of (lon, lat) coords, where N > 2. 
  (e.g. np.array([[-140,70], [-150, 75], [-145, 72]]))
  * There is no need to 'close' the polygon by repeating the final coordinate
  
- ellipsoid: named ellipsoid used to create polygon 
  (default: "WGS84")
  
- quiet: bool True/False, whether or not to suppress print statements as function runs 
  (default: True)

OUTPUT:
- poly: Shapely polygon

DEPENDENCIES:
import numpy as np, numpy.ma as ma
from pyproj import Geod
from shapely import wkt

Latest recorded update:
06-23-2023
    """
    
    
    # check that perim_coords was given with correct shape
    #*****************************************************
    assertation_print = f"perim_coords should have shape (N x 2), where N>2. "
    # shoul be two dimensional array
    assert len(np.shape(perim_coords))==2,assertation_print+f"Expected 2-d array, got {len(np.shape(perim_coords))}-dimensional array."
    # should be at least 3 points in array to make polygon
    assert np.shape(perim_coords)[0]>=3, assertation_print+f"Got ({np.shape(perim_coords)[0]} x {np.shape(perim_coords)[1]})"
    #*****************************************************
    

    # make polygon
    #-------------
    # specify a named ellipsoid
    geod = Geod(ellps=ellipsoid)
    
    # check whether polygon is already closed with given coordinates
    if perim_coords[0][0] == perim_coords[-1][0] and perim_coords[0][1] == perim_coords[-1][1]:
        perim_coords = perim_coords[0:-1]
    
    # run through coordinates and add to polygon
    WKT_STR = 'POLYGON (('
    for ii in range(len(perim_coords)):
        WKT_STR=WKT_STR+f'{str(perim_coords[ii,:][0])} {str(perim_coords[ii,:][1])},'
    # close polygon by repeating first coordinate
    WKT_STR=WKT_STR+f'{str(perim_coords[0,:][0])} {str(perim_coords[0,:][1])}))'

    # generate polygon
    poly = wkt.loads(WKT_STR)
        
    if not quiet:
        print(f'--> created: {WKT_STR}')
        area = abs(geod.geometry_area_perimeter(poly)[0])
        print('--> Geodesic area: {:12.1f} sq. m'.format(area))
        print('-->                {:12.1f} sq. km'.format(area/1e6))
    
    return poly



#///////////////////////////////////
#  make_smooth_geographic_box   ///
#/////////////////////////////////
#---------------------------------------------------------------------
# create smooth geographic box that smoothly follows parallels and 
# meridians between input coordinates
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np, numpy.ma as ma
#---------------------------------------------------------------------

def make_smooth_geographic_box(BOUND_LAT=[70.5,78], BOUND_LON=[187,220], num_points = 10):
    
    """Create smooth geographic box that smoothly follows parallels and 
meridians between input corner coordinates.

INPUT: 
- BOUND_LAT: [South, North] boundaries of box (default: [70.5,78])
- BOUND_LON: [West, East] boundaries of box (default: [187,220])
- num_points: extra points to add between longitude bounds to smooth curves (default: 10)

OUTPUT:
- box_lons: array of longitude values
- box_lats: array of latitude values

DEPENDENCIES:
import numpy as np, numpy.ma as ma

Latest recorded update:
03-28-2022
    """
    
    
    # check that BOUND_LAT, BOUND_LON were given with correct shape
    #*****************************************************
    assertation_print = f"BOUND_LAT and BOUND_LON should have shape (1 x 2). "
    assert len(BOUND_LAT)==2,assertation_print+f"Got BOUND_LAT with length {len(BOUND_LAT)}."
    assert len(BOUND_LON)==2,assertation_print+f"Got BOUND_LON with length {len(BOUND_LON)}."
    #*****************************************************

    
    # make smooth geographic box
    #---------------------------
    # initiate straight line along lon0 boundaries from lat0 to lat1
    box_lons = np.array([BOUND_LON[0],BOUND_LON[0]])
    box_lats = np.array([BOUND_LAT[0],BOUND_LAT[1]])
    
    # add num_points steps between lon0 and lon1 along lat1
    for ii in range(1,num_points+1):
        box_lons = np.append(box_lons, BOUND_LON[0]+(ii*(BOUND_LON[1]-BOUND_LON[0])/(num_points+1)))
        box_lats = np.append(box_lats, BOUND_LAT[1])
        
    # add extra straight line along lon1 boundaries from lat1 to lat0
    box_lons = np.append(box_lons,BOUND_LON[1])
    box_lons = np.append(box_lons,BOUND_LON[1])
    box_lats = np.append(box_lats,BOUND_LAT[1])
    box_lats = np.append(box_lats,BOUND_LAT[0])
    
    # add num_points steps between lon1 and lon0 along lat0
    for ii in range(1,num_points+1):
        box_lons = np.append(box_lons, BOUND_LON[1]+(ii*(BOUND_LON[0]-BOUND_LON[1])/(num_points+1)))
        box_lats = np.append(box_lats,BOUND_LAT[0])
        
    # close box at initial point
    box_lons = np.append(box_lons,BOUND_LON[0])
    box_lats = np.append(box_lats,BOUND_LAT[0])
    
    return box_lons, box_lats



#///////////////////////
#  make_linestring  ///
#/////////////////////
#---------------------------------------------------------------------
# create shapely linsetring from array of coordinates
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
from shapely import wkt
#---------------------------------------------------------------------

def make_linestring(coords):
    
    """Create shapely linsetring from (Nx2) array of coordinates.

INPUT: 
- coords: Nx2 array of (lon, lat) coords, where N > 2. 
  (e.g. np.array([[-140,70], [-150, 75], [-145, 72]]))

OUTPUT:
- LINESTRING: shapely linestring

DEPENDENCIES:
import numpy as np
from shapely import wkt

Latest recorded update:
03-28-2022
    """
    
    # check that coords was given with correct shape
    #*****************************************************
    # shoul be two dimensional array
    assert len(np.shape(coords))==2,f"coords should have shape (N x 2). Expected 2-d array, got {len(np.shape(coords))}-dimensional array."
    #*****************************************************
    
    # make linstring
    #---------------
    # run through coordinates and add to polygon
    WKT_STR = f'LINESTRING ({str(coords[0,:][0])} {str(coords[0,:][1])}'
    for ii in range(len(coords)):
        WKT_STR=WKT_STR+f', {str(coords[ii,:][0])} {str(coords[ii,:][1])}'
    WKT_STR=WKT_STR+')'
    
    # generate LINESTRING
    LINESTRING = wkt.loads(WKT_STR)
    
    return LINESTRING


#///////////////////////
#  regrid_geo_data  ///
#/////////////////////
#---------------------------------------------------------------------
# regrid data given initial data grid and desired geo grid
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
from scipy.interpolate import griddata
#---------------------------------------------------------------------

def regrid_geo_data(latgrid_initial,longrid_initial,datagrid_initial,latgrid_final,longrid_final, 
                   check_nan_grid = 'None', regrid_method = 'linear'):
    """Regrid data given initial data grid and desired geo grid.

INPUT: 
- latgrid_initial: (M x N) array of latitudes from initial geo grid
- longrid_initial: (M x N) array of longitudes from initial geo grid
- datagrid_initial: (M x N) array of data values from initial geo grid
- latgrid_final: (m x n) array of latitudes from desired geo grid
- longrid_final: (m x n) array of longitudes from desired geo grid
- check_nan_grid: if desired, supply grid that may contain nans that matches shape of input grids.
  If none supplied, will default to checking datagrid_initial for nans and eliminating them
- regrid_method: method for regridding data (defualt: 'linear')
                   

OUTPUT:
- datagrid_final: (m x n) array of regridded data on desired geo grid

DEPENDENCIES:
import numpy as np
from scipy.interpolate import griddata

Latest recorded update:
05-03-2022
    """
    
    
    # if initial data grid is 2D, concatenate along axis to 1D array
    if len(np.shape(datagrid_initial)):
        reshape_data = np.concatenate(datagrid_initial,axis=0)
        reshape_lat_initial = np.concatenate(latgrid_initial,axis=0)
        reshape_lon_initial = np.concatenate(longrid_initial,axis=0)
        # check for nans to eliminate from calcs
        if str(check_nan_grid) != 'None':
            print('Using nan_grid to eliminate nans from calculations.')
            NO_NANS = (np.isnan(np.concatenate(check_nan_grid,axis=0)) == False)
        else:
            print('Checking for nans in datagrid_initial to eliminate nans from calculations.')
            NO_NANS = (np.isnan(reshape_data) == False)
        
    # if already 1D, do nothing
    elif len(np.shape(datagrid_initial)):
        reshape_data = datagrid_initial
        reshape_lat_initial = latgrid_initial
        reshape_lon_initial = longrid_initial
        # check for nans to eliminate from calcs
        if str(check_nan_grid) != 'None':
            print('Using nan_grid to eliminate nans from calculations.')
            NO_NANS = (np.isnan(check_nan_grid) == False)
        else:
            print('Checking for nans in datagrid_initial to eliminate nans from calculations.')
            NO_NANS = (np.isnan(reshape_data) == False)
    else:
        print('Unfamiliar array dimensions of input grids. Must be 1D or 2D.')
    
    # eliminate nans from initial grids if they are present
    reshape_data = reshape_data[NO_NANS]
    reshape_lat_initial = reshape_lat_initial[NO_NANS]
    reshape_lon_initial = reshape_lon_initial[NO_NANS]
    
    # generate input coords to reshape
    INPUT_COORDS = np.stack((reshape_lon_initial,reshape_lat_initial),axis=1)
        
    # regrid 
    datagrid_final = griddata(INPUT_COORDS, reshape_data, (longrid_final, latgrid_final), method=regrid_method)

    return datagrid_final





        
#//////////////////////////////
#  within_polygon_indices  ///
#////////////////////////////
#---------------------------------------------------------------------
# Find indices of geo grid that fall within polygon.
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
#---------------------------------------------------------------------

def within_polygon_indices(polygon, lat_grid, lon_grid, quiet = True):
    
    """Find indices of geo grid that fall within polygon. Useful when repeatedly calculating mean of data values in same polygon as finding within-polygon indices is time-intensive.
    
INPUT:
- polygon: Shapely polygon
- lat_grid: MxN latitude grid
- lon_grid: MxN longitude grid
- quiet: bool True/False, whether or not to suppress print statements as function runs
    (default: True)

OUTPUT:
- poly_indices: N x 2 array of within-polygon indices for lat_grid/lon_grid. N [ii,jj] pairs to access lat_grid[ii,jj] and lon_grid[ii,jj] points.

DEPENDENCIES:
import numpy as np
import numpy.ma as ma
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

Latest recorded update:
06-23-2023
    """
    
    # run through coordinates in geographic grids and determine whether they're in polygon
    poly_Coords = []
    for ii in range(np.shape(lat_grid)[0]):
        for jj in range(np.shape(lat_grid)[1]):
            current_lon = lon_grid[ii,jj]
#             if current_lon>180:
#                 current_lon -= 360
            point = Point(current_lon,lat_grid[ii,jj])
            if polygon.contains(point):
                poly_Coords.append(ii)
                poly_Coords.append(jj)
    poly_coords = np.reshape(np.array(poly_Coords),(int(len(poly_Coords)/2),2))

    # find longitude/latitude within each box, and variables of interest
    poly_lon = np.array([])
    poly_lat = np.array([])
    poly_variable = np.array([])
    for ii in range(len(poly_coords)):
        poly_lon = np.append(poly_lon,lon_grid[poly_coords[ii][0],poly_coords[ii][1]])
        poly_lat = np.append(poly_lat,lat_grid[poly_coords[ii][0],poly_coords[ii][1]])
    
    if not quiet:
        print()
        print('Using shapely to find points in polygon')
        print(f'{len(poly_lon)} points found')
    
    # save polygon corner coordinates (ensuring lons are (0,360) to remove discontinuities in Alaskan Arctic)
    perimeter_lons = []
    perimeter_lats = []
    for ii in range(len(polygon.exterior.coords)):
        if polygon.exterior.coords[ii][0]<0:
            perimeter_lons.append(polygon.exterior.coords[ii][0]+360)
        else:
            perimeter_lons.append(polygon.exterior.coords[ii][0])
        perimeter_lats.append(polygon.exterior.coords[ii][1])
    
    # do a manual check on within-polygon points
    point_coords_to_remove = []
    
    for ii in range(len(poly_lon)):
        current_lon = poly_lon[ii]
        current_lat = poly_lat[ii]
        if current_lon<0:
            current_lon += 360
        if current_lon<np.min(perimeter_lons)-5 or current_lon>np.max(perimeter_lons)+5 or current_lat>np.max(perimeter_lats)+1 or current_lat<np.min(perimeter_lats)-1:
            point_coords_to_remove.append(ii)
    
    if not quiet:
        print()
        print('Checking shapely within-polygon points')
        print(f'--> Removing {len(point_coords_to_remove)} points that do not fall in/around polygon')
        print('....')
    
    # do a manual check on within-polygon points
    poly_lon_corrected = np.array([])
    poly_lat_corrected = np.array([])
    poly_coords_corrected = np.array([])

    for II in range(len(poly_lon)):
        if II not in point_coords_to_remove:
            poly_lon_corrected = np.append(poly_lon_corrected, poly_lon[II])
            poly_lat_corrected = np.append(poly_lat_corrected, poly_lat[II])
            poly_coords_corrected = np.append(poly_coords_corrected, poly_coords[II])

    # reshape indices of coords within polygon
    reshaped_coords = np.reshape(poly_coords_corrected, (int(len(poly_coords_corrected)/2),2))
            
    if not quiet:
        print(f'{len(poly_lon_corrected)} points now found within polygon')
        print()

    poly_indices = reshaped_coords.astype(int)
    
    return poly_indices




#/////////////////////////////////
#  polygon_mean_from_indices  ///
#///////////////////////////////
#---------------------------------------------------------------------
# Find mean value of data within polygon from indices of data_grid
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------

def polygon_mean_from_indices(polygon_indices, data_grid, lat_grid, weight_lats=False, quiet = True):
    
    """Find mean value of data within polygon from indices of data_grid (and associated lat_grid/lon_grid) that fall within polygon. Apply latitude-weighting from lat_grid if geo grids do not have equidistant spacing.
    
INPUT:
- polygon_indices: N x 2 array of within-polygon indices for lat_grid/lon_grid. N [ii,jj] pairs to access lat_grid[ii,jj] and lon_grid[ii,jj] points.
- data_grid: MxN variable grid
- lat_grid: MxN latitude grid
- weight_lats: bool True/False, whether or not to weight mean by latitudes
        (for use in data gridded in non equal-area projections)
        (default: False)
- quiet: bool True/False, whether or not to suppress print statements as function runs
    (default: True)

OUTPUT:
- poly_variable_mean: mean value of data_grid within given polygon

DEPENDENCIES:
import numpy as np
import numpy.ma as ma

Latest recorded update:
06-23-2023
    """
    
    polygon_data = np.array([])
    polygon_lats = np.array([])
    
    # run through and grab data
    for coordinate in polygon_indices:
        polygon_data = np.append(polygon_data, data_grid[coordinate[0],coordinate[1]])
        polygon_lats = np.append(polygon_lats, lat_grid[coordinate[0],coordinate[1]])
        
    if not quiet:
        print('Calculating variable mean within polygon')
        
    if weight_lats==False:
        if not quiet:
            print('--> Latitude weighting NOT applied')
        poly_variable_mean = np.nanmean(polygon_data)
    else:
        if not quiet:
            print('--> Latitude weighting applied')
        lat_weights=np.cos(np.pi*polygon_lats/180)
        poly_variable_mean = np.sum(lat_weights*polygon_data)/np.sum(lat_weights)
        
    if not quiet:
        print('....')
        print(f'variable mean: {poly_variable_mean:.2f}')
        
    return poly_variable_mean
        