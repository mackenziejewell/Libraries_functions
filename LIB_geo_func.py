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
#//////////////////////////////
#  distance_weighted_mean  ///
#////////////////////////////
#---------------------------------------------------------------------
# Find inverse distance - weighted mean value of data at provided lat/lon coordinate from gridded data.
#---------------------------------------------------------------------
#///////////////////////////
#  make_geodesic_paths  ///
#/////////////////////////
#---------------------------------------------------------------------
# Create geodesic paths
#---------------------------------------------------------------------
#//////////////////////
#  lonlat_to_proj  ///
#////////////////////
#---------------------------------------------------------------------
# Convert lat/lon coordinates to projected coordinates using cartopy.
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
        
    
    
#//////////////////////////////
#  distance_weighted_mean  ///
#////////////////////////////
#---------------------------------------------------------------------
# FFind inverse distance - weighted mean value of data at provided lat/lon coordinate from gridded data.
#---------------------------------------------------------------------
# DEPENDENCIES:
from pyproj import Geod
import numpy as np
from metpy.units import units
#---------------------------------------------------------------------

def distance_weighted_mean(point_lon = [], point_lat = [], grid_lons = [], grid_lats = [], 
                       grid_data = [], mask = [],  ellps = 'WGS84',
                       lat_buffer = 1, lon_buffer = 5,
                       return_vars = ['neighbor_lons', 'neighbor_lats', 'neighbor_data', 'neighbor_weights', 'weighted_mean'],
                       max_dist = 50*units('km')):
    
    """Find inverse distance - weighted mean value of data at provided lat/lon coordinate from gridded data. Option to simply find and return data from nearest gridded point.
    
INPUT:
- point_lon: longitude to find nearest on grid (range: 0,360)
- point_lat: latitude to find nearest on grid
- grid_lons: M x N array of longitude values for geospatial grid (range: 0,360)
- grid_lats: M x N array of latitude values for geospatial grid
- grid_data: M x N array of data values corresponding to geospatial grid
- mask: M x N array filled with bools, True for masked regions to not be considered, and False for
        gridded lons/lats to include in distance calculations. 
        If no mask desired, set to [] (default)
- ellps: named ellipsoid used to create polygon (default: "WGS84")
- lat_buffer: maximum allowed difference in degrees between latitudes of input vs gridded coordinates.
    Check this first, if within lat_buffer degrees, then calculate geodesic distance. This is simply a time-saver. 
    (default: 1 degrees)
- lon_buffer: maximum allowed difference in degrees between longitude of input vs gridded coordinates.
    Check this first, if within lon_buffer degrees, then calculate geodesic distance. This is simply a time-saver. 
    (default: 5 degrees)
    
- max_dist: Pint quantity (with Metpy units), maximum allowed distance between provided and nearest gridded point
    (default: 50*units('km') for 50 km). If nearest gridded point is further than this, return empty values 
    for nearest points.


OUTPUT:
List of any or all of variables specified in return_vars:
- nearest_lon: nearest gridded lon to provided point, empty list if none within max_dist
- nearest_lat: nearest gridded lat to provided point, empty list if none within max_dist
- nearest_dist: distance of nearest gridded point to provided point, flag value 999.0 km if none within max_dist
                if nearest_dist == 0.0 (direct match to gridded point), all output variables with 
                'neighbor' and 'weight' in name will be filled with single values from corresponding matching point.
- nearest_data: data of nearest gridded point to provided point, empty list if none within max_dist
- all_dist: M x N array of gridded distances (with metpy units) from provided point, 
            filled with value representing 999 km where data is masked, points outside lat/lon buffers,
            and for points farther than max_dist away
- neighbor_lons: 1 X D array of neighboring longitudes within max_dist of provided point
                (single value array if direct spatial match)
- neighbor_lats: 1 X D array of neighboring latitudes within max_dist of provided point
                (single value array if direct spatial match)
- neighbor_dist: 1 X D array of distances (with metpy units) between provided point and neighboring points within max_dist
                (single value array if direct spatial match)
- neighbor_data: 1 X D array of data values for gridded points within max_dist of provided point
                (single value array if direct spatial match)
- neighbor_weights: 1 X D array of inverse distance weights (with metpy units) for gridded points within max_dist of provided point
                (single value array of unity, unitless if direct spatial match)                
- weighted_mean: inverse distance-weighted average value of data across gridded points within max_dist of provided point
                (if direct spatial match, directly return data of corresponding grid point)  
        
        
DEPENDENCIES:
from pyproj import Geod
import numpy as np
from metpy.units import units

Latest recorded update:
07-05-2023
    """
    
    
    # check for input errors
    if point_lon < 0:
        raise ValueError(f'point_lon should be greater than or equal to zero (range: 0-360), not {point_lon}')
    if (grid_lons < 0).any():
        raise ValueError(f'grid_lons should be greater than or equal to zero (range: 0-360). {np.sum(grid_lons < 0)} values found below 0.')
    if str(type(max_dist)) !=  "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>":
        print(str(type(max_dist)))
        raise TypeError(f"max_dist should include metpy units. If units are km: set as {max_dist}*units('km')")
    assert len(return_vars) > 0, 'return_vars list is empty. Must have length >=1'
    
    # create Geod from ellipsoid name
    g = Geod(ellps=ellps)

    # use geodesics to find nearest grid coordinate to lead coordinate
    nearest_dist = max_dist
    
    nearest_lat = []
    nearest_lon = []
    nearest_data = []
    
    # if no mask provided, include all values
    if len(mask) == 0:
        mask = np.full(grid_lons.shape, False)
        
    # if no data provided, fill with nans
    if len(grid_data) == 0:
        grid_data = np.full(grid_lons.shape, np.nan)
        
    # create fill vall of 999 km in whatever units provided, then strip units 
    # create array to save distances from provided points
    fill_val = ( ( 999.0 * units('km') ).to(max_dist.units) ).magnitude
    all_dist = np.full_like(grid_lons, fill_val)
    
    # run through all gridded points to look for closest point
    for II in range(np.shape(grid_lons)[0]):   
        for JJ in range(np.shape(grid_lons)[1]):

            # only calculate distance for not masked values
            if mask[II,JJ] == False:

                # if grid location within buffer degrees of point
                # use pyproj to calculate actual distance
                if np.abs(point_lat-grid_lats[II,JJ])<lat_buffer and np.abs(point_lon-grid_lons[II,JJ])<lon_buffer:

                    # calculate distance between grid point and provided point
                    # pyproj runs slightly faster than geopy's geodesic            
                    az12,az21,distance = g.inv(point_lon, point_lat, grid_lons[II,JJ], grid_lats[II,JJ])
                    # convert distance to provided units, then strip units to save to array
                    distance = (distance*units('m')).to(max_dist.units)
                    all_dist[II, JJ] = distance.magnitude

                    # if closer than nearest_dist, replace nearest coords and dist
                    if distance < nearest_dist:
                        nearest_dist = distance
                        nearest_lat = grid_lats[II,JJ]
                        nearest_lon = grid_lons[II,JJ]
                        nearest_data = grid_data[II,JJ]

    # if no value found within max_dist of provided point, notify and flag nearest_dist
    if nearest_dist == max_dist:
        print(f'Could not find grid coordinate within {max_dist} of {point_lat, point_lon}')
        nearest_dist = ( 999.0 * units('km') ).to(max_dist.units)
        
    # store all output in dictionary
    vars_dict = {}
    vars_dict['nearest_lon'] = nearest_lon
    vars_dict['nearest_lat'] = nearest_lat
    vars_dict['nearest_dist'] = nearest_dist
    vars_dict['nearest_data'] = nearest_data
    
    # replace points further than max_dist away with fill value
    # save to dict with units
    all_dist[all_dist > max_dist.magnitude] = fill_val
    vars_dict['all_dist'] = all_dist * max_dist.units
    
    
    # if there is a direct match with nearest gridded point and value, 
    # save nearest points found above as single-value array
    if nearest_dist.magnitude == 0.0:
        vars_dict['neighbor_lons'] = np.array([nearest_lon])
        vars_dict['neighbor_lats'] = np.array([nearest_lat])
        vars_dict['neighbor_dist'] = np.array([nearest_dist.magnitude])*nearest_dist.units
        vars_dict['neighbor_data'] = np.array([nearest_data])
        vars_dict['neighbor_weights'] = np.array([1])
        vars_dict['weighted_mean'] = nearest_data
        
    # if multiple values within max_dist of nearest point, 
    # calculate inverse distance - weighted mean of data 
    else:  
        
        # find indices of gridded points within max_dist range of point
        iii, jjj = np.where(all_dist < fill_val)

        # save points to 1-d arrays
        neighbor_lons = grid_lons[iii, jjj]
        neighbor_lats = grid_lats[iii, jjj]
        neighbor_dist = all_dist[iii, jjj]
        neighbor_data = grid_data[iii, jjj]

        # create inverse distance weights
        neighbor_weights = (1/neighbor_dist)
        
        # calculate weighted mean data value
        weighted_mean = np.sum(neighbor_data  * neighbor_weights) / np.sum(neighbor_weights)
        
        vars_dict['neighbor_lons'] = neighbor_lons
        vars_dict['neighbor_lats'] = neighbor_lats
        vars_dict['neighbor_dist'] = neighbor_dist * max_dist.units
        vars_dict['neighbor_data'] = neighbor_data
        vars_dict['neighbor_weights'] = neighbor_weights * (1/max_dist.units)
        vars_dict['weighted_mean'] = weighted_mean
        
    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    if len(return_data) == 1:
        return_data = return_data[0]
    return return_data



#///////////////////////////
#  make_geodesic_paths  ///
#/////////////////////////
#---------------------------------------------------------------------
# Create geodesic paths
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
from metpy.units import units
from pyproj import Geod
#---------------------------------------------------------------------

def make_geodesic_paths(method = 'constant_azimuth',
                        c1 = (200, 70), c2 = (200, 71), N = 10, 
                        azimuth = 0 * units('degree'), distance = 25*units('km'), 
                        g = Geod(ellps='WGS84'), quiet = True):
    
    """Create equally-spaced paths along constant geodesic between specified points or with constant azimuth.

INPUT: 
- method: method to create the path:
    if 'constant_azimuth': use c1, azimuth, and N to create path with N-1 new points from c1.
                            At each step, create next step from local azimuth orientation. 
                            Terminal point cannot be specified. (uses pyproj g.fwd)
    if 'geodesic_npts': use c1, c2, and N to create path with N-2 new points between c1 and c2.
                        Local azimuth changes at each point if lon/lat changes between c1/c2. 
                        Points are evenly-spaced, but spacing size and local azimuths cannot 
                        be specified. (uses pyproj g.npts)
- c1: tuple of (lon, lat) of starting point.
- c2: tuple of (lon, lat) of terminus point.
- N: int, total number of points in path. Must be 2 or greater.
- azimuth: angle in degrees or radians from Northward (+ CW)
           must include metpy units as degrees or radians (default 0 * units('degree'))
- g: geod to use for pyproj geodesic calculations (default: Geod(ellps='WGS84'))
- quiet: bool, whether or not to suppress prints (default: True)

OUTPUT:
- path: (N x 2) array of lon (axis 0) and lat (axis 1) coords. 

DEPENDENCIES:
import numpy as np
from metpy.units import units
from pyproj import Geod

Latest recorded update:
07-06-2023
    """
    
    if str(type(distance)) !=  "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>":
        print(str(type(distance)))
        raise TypeError(f"distance should include metpy units. If units are km: set as {distance}*units('km')")
    if str(type(azimuth)) !=  "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>":
        print(str(type(azimuth)))
        raise TypeError(f'azimuth should include metpy units "degree" or "radian". If units are degree: set as {azimuth}*units("degree")')
    if N < 2:
        raise ValueError('N must be 2 or greater.')
        
    # convert distance to meters and azimuth to degrees
    distance_m = ( distance.to('m') ).magnitude
    azimuth_d = ( azimuth.to('degree') ).magnitude
    
    # unpack coords
    lon1, lat1 = c1
    lon2, lat2 = c2
    
    if not quiet:
        print(f'>> create path with {method} method, N = {N} points') 
        print(f'>> start point: {c1}')
        
        if str(method) == 'geodesic_npts':
            print(f'>> terminus point: {c2}')
            print(f'>> add {N-2} points in between')
        elif str(method) == 'constant_azimuth':
            print(f'>> add {N-1} points, maintaining azimuth = {azimuth}')
        
    
    # create path including initial / final points
    #----------------------------------------------
    # initial point
    path_lon = np.array([lon1])
    path_lat = np.array([lat1])
    
    
    # middle and final points, depending on method
    
    # method 1
    #----------------------
    if str(method) == 'geodesic_npts':
        # create path with N equally spaced points between start and terminus
        lonlats = g.npts( lon1,lat1, lon2,lat2, N-2 )
        # intermediate points
        for ll in range(len(lonlats)):
            path_lon = np.append(path_lon, lonlats[ll][0])
            path_lat = np.append(path_lat, lonlats[ll][1])
        # final point
        path_lon = np.append(path_lon, lon2)
        path_lat = np.append(path_lat, lat2)
    
    # method 2
    #----------------------
    elif str(method) == 'constant_azimuth':
        for ii in range(N-1):
            if ii == 0:
                lon2, lat2 = lon1, lat1
            lon2, lat2, az = g.fwd(lon2, lat2, azimuth_d, distance_m)
            path_lon = np.append(path_lon, lon2)
            path_lat = np.append(path_lat, lat2)
            
    # convert to (0, 360) range
    path_lon[path_lon<0]+=360

    # stack 
    path = np.stack((path_lon, path_lat), axis=1)
    
    return path




#//////////////////////
#  lonlat_to_proj  ///
#////////////////////
#---------------------------------------------------------------------
# Convert lat/lon coordinates to projected coordinates using cartopy.
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
#---------------------------------------------------------------------

def lonlat_to_proj(lonlats = np.array([[200,70],[200,71]]), proj_crs=[], quiet = True):
    
    """Convert lat/lon coordinates to projected coordinates using cartopy.
    
INPUT:
- lonlats: N x 2 array of (lon, lat) coordinates. Lons are axis 0, lats are axis 1.
    default: np.array([[200,70],[200,71]])
- proj_crs: cartopy projection 
    (e.g. ccrs.LambertAzimuthalEqualArea(central_longitude=0, central_latitude=90,
    globe=ccrs.Globe(semimajor_axis = 6371228, semiminor_axis = 6371228)))
- quiet: bool True/False, whether or not to suppress print statements as function runs 
  (default: True)
  
OUTPUT:
- xx_yy_proj: N x 2 array of projected (xx, yy) coordinates

DEPENDENCIES:
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs

Latest recorded update:
07-06-2023
    
    """
    
    if type(proj_crs) == list:
        raise TypeError(f'proj_crs should be cartopy.crs projection, not {type(proj_crs)}')
    
    if not quiet:
            print(f'(lon, lat) --> (xx, yy):\n-----------------------')
                  
    # transform point by point and save to xx_yy_projs
    xx_yy_projs = np.array([])
    for lonlat in lonlats:
        # make the transformation
        xx_yy = proj_crs.transform_point(*lonlat, src_crs=ccrs.PlateCarree())
        xx_yy_projs = np.append(xx_yy_projs, xx_yy)
        if not quiet:
            print(f'({lonlat[0]:.2f}, {lonlat[1]:.2f}) --> ({xx_yy[0]:.2f}, {xx_yy[1]:.2f})')
        
    # reshape to be (N X 2)
    xx_yy_projs = np.reshape(xx_yy_projs, np.shape(lonlats))
    
    return xx_yy_projs
    