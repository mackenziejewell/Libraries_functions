#///////////////////////
#  gate_drift_comp  ///
#/////////////////////
#---------------------------------------------------------------------
# Find ice drift components relative to flux gates.
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
from metpy.units import units
from LIB_geo_func import (distance_weighted_mean)
#---------------------------------------------------------------------

def gate_drift_comp(gate = [], azimuths = [], u_grid = [], v_grid = [], lat_grid = [], lon_grid = [],
                    max_dist = 50*units('km'), 
                    return_vars = ['gate_u', 'gate_v', 'perp_vecs', 'para_vecs', 'gate_perp', 'gate_para'],
                   quiet = False):
    
    """Find ice drift components relative to flux gates.
    
INPUT:
- gate: (N x 2) array of lon (axis 0) and lat (axis 1) coords. 
- azimuth: (N x 1) array of gate angles in degrees or radians from Northward (+ CW)
           must include metpy units as degrees or radians (default 0 * units('degree'))
- u_grid: (M x L) grid of u (eastward) ice drift components
- v_grid: (M x L) grid of v (northward) ice drift components
- lat_grid: (M x L) grid of latitudes corresponding to u_grid, v_grid
- lon_grid: (M x L) grid of longitude corresponding to u_grid, v_grid
- max_dist: Pint quantity (with Metpy units), maximum allowed distance between provided and nearest gridded point
    (default: 50*units('km') for 50 km). If nearest gridded point is further than this, return empty values 
    for nearest points.
- return_vars: list of variable names to return
 - quiet: bool, whether or not to suppress prints (default: False)
 
OUTPUT:
List of any or all of variables specified in return_vars:
- gate_u: N x 1 array of eastward drift components at gate points
- gate_v: N x 1 array of northward drift components at gate points
- perp_vecs: (N x 2) array of gate-perpendicular (90 CW from gate-parallel) vectors (length 1)
            (axis 0, axis 1) as (perp_i, perp_j) 
            with perp_i eastward, perp_j northward. 90 deg CCW from gate-parallel.
- para_vecs: (N x 2) array of gate-parallel (along-azimuth) vectors (unit length)
            (axis 0, axis 1) as (para_i, para_j) 
            with para_i eastward, para_j northward. 90 deg CW from gate-perpendicular.
- gate_perp: gate-normal ice drift components
- gate_para: gate-parallel ice drift components

DEPENDENCIES:
import numpy as np
from metpy.units import units
sys.path.append('../Libraries_functions/')
from LIB_geo_func import (distance_weighted_mean)

# homemade function from: GitHub/Libraries_functions/LIB_geo_func.py
distance_weighted_mean

Latest recorded update:
07-06-2023
    
    """
    
    
    if str(type(azimuths[0])) !=  "<class 'pint.quantity.build_quantity_class.<locals>.Quantity'>":
        raise TypeError(f"azimuths should include metpy units. If units are degrees: set as azimuths*units('degree')")
    
    # empty arrays to store data
    gate_u = np.array([])
    gate_v = np.array([])
    
    perp_vecs = np.array([])
    para_vecs = np.array([])
    
    gate_perpendicular = np.array([])
    gate_parallel = np.array([])
    
    gate_norm_x = np.array([])
    gate_norm_y = np.array([])
    
    # find drift at all gate points
    #------------------------------
    for ii in range(len(gate)):
    
        # use distance weighted mean to approx local drift
        # LATE CHANGE TO INTERPOLATION
        #--------------------------------------------------
        # U component
        u_mean = distance_weighted_mean(point_lon = gate[ii, 0], point_lat = gate[ii, 1], 
                                        grid_lons = lon_grid, grid_lats = lat_grid, 
                                        grid_data = u_grid, mask = np.isnan(u_grid),  
                                        return_vars = ['weighted_mean'],
                                        max_dist = max_dist, quiet = True)
        # V component
        v_mean = distance_weighted_mean(point_lon = gate[ii, 0], point_lat = gate[ii, 1], 
                                        grid_lons = lon_grid, grid_lats = lat_grid, 
                                        grid_data = v_grid, mask = np.isnan(v_grid),  
                                        return_vars = ['weighted_mean'],
                                        max_dist = max_dist, quiet = quiet)

        # save eastward, northward drift components
        gate_u = np.append(gate_u, u_mean)
        gate_v = np.append(gate_v, v_mean)

        # create unit vectors for gate
        #-----------------------------
#         # gate-perpendicular
#         gate_perp_i = -np.cos(azimuths[ii].to('radian').magnitude)
#         gate_perp_j = np.sin(azimuths[ii].to('radian').magnitude)
#         perp_vecs = np.append(perp_vecs, (gate_perp_i, gate_perp_j))
        
#         # gate-parallel
#         # convert azimuth to beta coordinate
#         beta = (np.pi/2)*units('radian') - azimuths[ii].to('radian')
#         # ensure within [-180, 180] range
#         if beta <= -np.pi*units('radian'):
#             beta += 2*np.pi*units('radian')
#         gate_para_i = np.cos(beta.to('radian').magnitude)
#         gate_para_j = np.sin(beta.to('radian').magnitude)
#         para_vecs = np.append(para_vecs, (gate_para_i, gate_para_j))
        
    
        # gate-parallel
        gate_para_i = np.sin(azimuths[ii].to('radian').magnitude) # northward
        gate_para_j = np.cos(azimuths[ii].to('radian').magnitude) # eastward
        para_vecs = np.append(para_vecs, (gate_para_i, gate_para_j))
        
     
        # gate-perpendicular
        gate_perp_i = gate_para_j
        gate_perp_j = -gate_para_i
        perp_vecs = np.append(perp_vecs, (gate_perp_i, gate_perp_j))
        
        # find gate-perp and -para drift
        #-------------------------------
        # find gate-normal ice drift component
        gate_perp_drift = (u_mean * gate_perp_i) + (v_mean * gate_perp_j)
        gate_perpendicular = np.append(gate_perpendicular, gate_perp_drift)
        
        # find gate-parallel ice drift component
        gate_para_drift = (u_mean * gate_para_i) + (v_mean * gate_para_j)
        gate_parallel = np.append(gate_parallel, gate_para_drift)
        
    # store all output in dictionary
    vars_dict = {}
    
    vars_dict['gate_u'] = gate_u
    vars_dict['gate_v'] = gate_v
    
    perp_vecs = np.reshape(perp_vecs, np.shape(gate))
    para_vecs = np.reshape(para_vecs, np.shape(gate))
    vars_dict['perp_vecs'] = perp_vecs
    vars_dict['para_vecs'] = para_vecs
    
    vars_dict['gate_perp'] = gate_perpendicular
    vars_dict['gate_para'] = gate_parallel
    
    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    if len(return_data) == 1:
        return_data = return_data[0]
    return return_data
