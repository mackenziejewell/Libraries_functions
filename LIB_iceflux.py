

#////////////////////////
#  gate_interp_nans  ///
#//////////////////////
#---------------------------------------------------------------------
# Linearly interpolate across nan values in gate-orthogonal drift data
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
#---------------------------------------------------------------------
def gate_interp_nans(gate_perp_values):
    
    """Function to linearly interpolate across nans in gate-orthogonal drift data
    
INPUT:
- gate_perp_values: (N x 1) array of gate-orthogonal drift data.
    (should not have nans at endpoints! Use "gate_handle_nans" to ensure that first.)

OUTPUT:
- gate_perp: (N x 1) array of updated gate-orthogonal drift data, with interpolations across nan values


DEPENDENCIES:
import numpy as np

Latest recorded update:
04-02-2024
    
    """
    
    
    gate_perp = np.copy(gate_perp_values)
    
    assert gate_perp[0] != np.nan, "First value is np.nan! Should not have nans at endpoints."
    assert gate_perp[-1] != np.nan, "Last value is np.nan! Should not have nans at endpoints."
    
    # loop through all but endpoints
    for ii in range(1,len(gate_perp)-1):

        # if any nans are found, interpolate values from neighbors
        if np.isnan(gate_perp[ii]):
            n_e = ii

            # find non-nan neighbors
            #---------------------------
            # step backward to find nearest non-nan neighbore BEFORE error
            n_b = np.nan
            n_a = np.nan
            for jj in range(1,n_e)[::-1]:
                if not np.isnan(gate_perp[jj]):
                    n_b = jj
                    break
            # step forward to find nearest non-nan neighbore AFTER error
            for jj in range(n_e+1,len(gate_perp)):
                if not np.isnan(gate_perp[jj]):
                    n_a = jj
                    break
            if np.isnan(n_b) or np.isnan(n_a):
                print('Cant find any non-nan neighbors!')

            # interpolate value
            #---------------------------
            # dist between non-nans (b->a) and (b->e)
            dS = n_a - n_b
            ds = n_e - n_b

            # difference between values (b->a) and (b->e)
            dV = gate_perp[n_a] - gate_perp[n_b]
            dv = (dV/dS)*ds

            # estimated values at (e)
            v_e = gate_perp[n_b] + dv

            # replace values
            gate_perp[n_e] = v_e

    return gate_perp




#////////////////////////
#  gate_handle_nans  ///
#//////////////////////
#---------------------------------------------------------------------
# Handle conditions/nans at individual gate points.
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
from metpy.units import units
#---------------------------------------------------------------------
def gate_handle_nans(gate_perp_values, set_zero = [], set_nan_zero = []):
    
    """Function to handle conditions/nans at individual gate points
    
INPUT:
- gate_perp_values: (N x 1) array of gate-orthogonal drift component data. 
- set_zero: indices of gate_perp_values to unconditionally set to zero
    (for use at gate endpoints ending on land)
- set_nan_zero: indices of gate_perp_values to set to zero IF they are np.nan
    (for use at gate endpoints ending on ice that is usually landfast)

OUTPUT:
- gate_perp: (N x 1) array of updated gate-orthogonal drift component data


DEPENDENCIES:
import numpy as np
from metpy.units import units

Latest recorded update:
04-02-2024
    
    """
    
    gate_perp = np.copy(gate_perp_values)
    
    # set to zero (should always be nan anyway)
    # use for points that are on land (gate endpoints)
    for ii in set_zero:
        gate_perp[ii] = 0 * units('cm/s')
    
    # set to zero if is a nan
    # use for points in LF ice. Likely always a nan (in which case will replace with 0 for stationary ice)
    # but allow values there in case it does ever move
    for ii in set_nan_zero:
        if np.isnan(gate_perp[ii]):
            gate_perp[ii] = 0 * units('cm/s')

    return gate_perp


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
