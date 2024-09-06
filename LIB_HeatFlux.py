

from metpy.units import units
import numpy as np
from metpy.calc import relative_humidity_from_dewpoint, mixing_ratio_from_relative_humidity, vapor_pressure
from scipy.interpolate import RegularGridInterpolator

# turbulent_OHF
# calc_friction_u
# QERA_atpoints
# calc_ssh
# calc_lwd
# calc_lwu
# seasonal_Cio
# calc_vars_from_T_S

# OLD FUNCTIONS
# calc_friction_v


def turbulent_OHF(dT, us, rho_o = 1023 * units('kg/m3'), 
                  Cp = 3980 * units('J') * units('kg')**(-1) * units('delta_degC')**(-1), 
                 St = 0.0057):

    """Calculate bulk turbulent heat flux. Reference McPhee 2008 and Zhong (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021GL096216)

INPUT: 
- dT: mixed layer temperature difference from freezing (units: 'delta_degC')
- us: surface friction velocity (units: 'm/s')
- rho_o: seawater density (default: 1023 * units('kg/m3'))
- Cp: seawater heat capacity (default: 3980 * units('J') * units('kg')**(-1) * units('delta_degC')**(-1))
- St: heat transfer coefficient / stanton number (default: St = 0.0057)

OUTPUT:
- F_H_bulk: estimated turbulent heat flux (units: 'W/meter^2')

DEPENDENCIES:
from metpy.units import units
import numpy as np

Latest recorded update:
09-05-2024
    """
    
    # bulk turbulent heat flux
    F_H_bulk = (rho_o * Cp * St * us * dT).to('W/meter^2')

    return F_H_bulk


def calc_friction_u(delta_u, z, z0 = 0.05 * units('m')):
    
    """Calculate friction velocity, based on McPhee 2008.

INPUT: 
- delta_u: magnitude difference between ice and ocean velocities (but retain units), with depth
- z0: roughness length scale
- z: ocean depth, negative from surface = 0

OUTPUT:
- us0: estimate surface friction velocity from ocean currents at depth

DEPENDENCIES:
from scipy.interpolate import RegularGridInterpolator
import numpy as np

Latest recorded update:
09-05-2024
    """
    
    K = 0.4 # van karman's constant
    
    # square root of effective drag coefficient 
    # with depth from effective drag length scale
    sq_cd = ( (1/K) * np.log(-z/z0) ) ** (-1)
    
    us0 = sq_cd * delta_u
    
    return us0




def QERA_atpoints(x_grid, y_grid, z_grid, x_pts, y_pts, method = 'linear'):
    
    """Interpolate ERA5 lat/long data to desired lat/lon coordinates.

INPUT: 
- x_grid: gridded longitudes
- y_grid: gridded latitudes
- y_grid: gridded data (no units!)
- x_pts: longitude of desired coordinates
- y_pts: latitude of desired coordinates
- method: interpolation method (default: 'linear')

OUTPUT:
- z_pts: intepolated data values

DEPENDENCIES:
from scipy.interpolate import RegularGridInterpolator
import numpy as np

Latest recorded update:
09-05-2024
    """
    
    #from scipy.interpolate import RegularGridInterpolator
    # transpose EAR5 data since its on backwards-seeming grid
    # convert all data to float64, says "no matching signature" if some are float32
    x = np.transpose(x_grid)[:,0].astype(np.float64)
    y = np.transpose(y_grid)[0,:].astype(np.float64)
    z = np.transpose(z_grid).astype(np.float64)

    # Create the interpolator
    interp = RegularGridInterpolator((x, y), z, method=method)

    # coordinates at which we want data
    points = np.array(list(zip(x_pts, y_pts)))
    z_pts = interp(points)
    
    return z_pts



def calc_ssh(Ta, Ts, p, U,
            Cpa = 1005 * units('J/(kg delta_degC)'),
            Cs = 0.00175):
    
    """Calculate turbulent air-surface sensible heat flux (+ upwards).
    Referencing Weeks: On Sea Ice text, page 198.

INPUT: 
- Ta: near-surface air temperature (include units: K or degC)
- Ts: surface temperature (include units: K or degC)
- p: surface pressure (include units: Pa)
- U: surface wind speed (include units: m/s)

constants with adjustable values:
- Cpa: specific heat of air (include units:  J kg-1 K-1)
- Cs: bulk transfer coefficient (assume 0.00175 from Weeks p.198)

OUTPUT:
- FS: surface sensible heat flux [W m-2] (+ is surface -> atmosphere)

DEPENDENCIES:
from metpy.units import units
import numpy as np

Latest recorded update:
09-03-2024
    """
    
    # specific gas constant for dry air
    Rspec = 287.0500676 * units('J/(kg delta_degC)').to('N m/(kg K)')
    
    # convert temperatures
    if Ta.units == units('degree_Celsius'):
        TA = Ta.to('K')
    elif Ta.units == units('K'):
        TA = Ta
    else:
        print(f'Unrecognized units for Ta: {Ta.units}')
        
    if Ts.units == units('degree_Celsius'):
        TS = Ts.to('K')
    elif Ts.units == units('K'):
        TS = Ts
    else:
        print(f'Unrecognized units for Ts: {Ts.units}')
        
    # calculate air density from surface pressure, temperature
    #---------------------------------------------------------
    # surface pressure
    rho = (p.to('N/m**2') / (Rspec * TA))

    # calculate surface heat flux
    FS = - rho * Cpa * Cs * U * (((TA - TS).magnitude)*units('delta_degC'))
    
    return FS.to('W/m**2')



def calc_lwu(Ts, Es = 0.98):
    
    """Calculate upward longwave radiation (+ upwards).
    Referencing Weeks: On Sea Ice text, page 197 and Pease 1987

INPUT: 
- Ts: surface temperature (include units: K or degC)

constants with adjustable values:
- Es : surface emissivity (default 0.98 from Pease 1987)

OUTPUT:
- LWU: upward longwave radiative flux [W m-2] (+ is surface -> atmosphere)

DEPENDENCIES:
from metpy.units import units
import numpy as np

Latest recorded update:
09-03-2024
    """
    
    # boltzmann constant
    SIG = 5.67 * 10**(-8) * units('W/(m**2 K**4)')
    
    # convert temperatures
    if Ts.units == units('degree_Celsius'):
        TS = Ts.to('K')
    elif Ts.units == units('K'):
        TS = Ts
    else:
        print(f'Unrecognized units for Ts: {Ts.units}')
        
    # calculate longwave upward flux
    LWU = Es * SIG * TS**4
    
    return LWU.to('W/m**2')


def calc_lwd(Ta, Td, p, C):
    
    """Calculate downward longwave radiation (+ upwards).
    Referencing Weeks: On Sea Ice text, page 197

INPUT: 
- Ta: near-surface air temperature (include units: K or degC)
- Td: near-surface dewpoint temperature (include units: K or degC)
- p: surface pressure (include units:Pa)
- e: near-surface partial vapor pressure of water (included units: Pa)
- C: cloud fraction

constants with adjustable values:
- Ea : atmospheric emissivity (default 0.95 from Pease 1987)

OUTPUT:
- LWD: downward longwave radiative flux [W m-2] (+ is surface -> atmosphere)

DEPENDENCIES:
from metpy.units import units
import numpy as np

Latest recorded update:
09-03-2024
    """

    # boltzmann constant
    SIG = 5.67 * 10**(-8) * units('W/(m**2 K**4)')
    
    # convert temperatures
    if Ta.units == units('degree_Celsius'):
        TA = Ta.to('K')
    elif Ta.units == units('K'):
        TA = Ta
    else:
        print(f'Unrecognized units for Ta: {Ta.units}')
        
    # calculate longwave upward flux (simple relation, using atmospheric emissivity 0.95 from Pease 1987)
    Ea = 0.95
    LWD_simple = - (Ea * SIG * TA**4).to('W/m**2')
    
    # calculate vapor pressure of water from Ta, Td, p
    rh = relative_humidity_from_dewpoint(Ta, Td).to('percent') # relative humidity
    mr = mixing_ratio_from_relative_humidity(p, Ta, rh).to('g/kg') # mixing ratio
    e = vapor_pressure(p, mr) # partial vapor pressure

    # LWD under clear-sky conditions 
    # use empirical estimate of Ea from Efimova 1961 (but ref.d from Weeks P. 197)
    Ea = (0.746 + 0.0066 * e.magnitude)
    LWD_clr = - (Ea * SIG * TA**4).to('W/m**2')
    
    # clear LWD to all LWD, from Maykut and Church 1973 (but ref.d from Weeks P. 198)
    LWD_all = LWD_clr * (1 + 0.22 * C ** (2.75))
    
    return LWD_simple, LWD_clr, LWD_all



# seasonal trend in ice-ocean drag constant from Brenner et al Figure 6
def seasonal_Cio(date):
    # estimate linear seasonal trend in Cio 
    # between Nov 1 and Feb 15

    # determine year corresponding to Jan/Feb in Nov - Feb seasonal cycle
    if date.month > 10:
        year = date.year + 1
    else:
        year = date.year

    # estimated from Brenner et al Figure 6
    # https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020JC016977
    date_f = datetime(year, 2, 15)
    Cf = 7 * 10**(-3)
    date_i = datetime(year-1, 11, 1)
    Ci = 3 * 10**(-3)

    # estimate current Cio value
    total_dt = (date_f - date_i).days
    current_dt = (date - date_i).days
    m = (Cf - Ci)/total_dt
    Cio = (m * current_dt) + Ci
    
    return Cio



import gsw
import gsw.freezing

def calc_vars_from_T_S(T, SP, depth, lon, lat, saturation_fraction = 1):
    
    # sea pressure ( i.e. absolute pressure - 10.1325 dbar )
    p = gsw.conversions.p_from_z(depth, lat, geo_strf_dyn_height=0, sea_surface_geopotential=0)

    # Absolute Salinity (g / kg)
    SA = gsw.conversions.SA_from_SP(SP, p, lon, lat)

    # Conservative Temperature (ITS-90) (degC)
    CT = gsw.conversions.CT_from_t(SA, T, p)

    # Saturation fraction of dissolved air in seawater. (0..1)
    # let's assume its 1?
#     saturation_fraction = 1
    
    # freezing point of water
    T_f = gsw.freezing.CT_freezing(SA, p, saturation_fraction)

    # potential density
    sigma0 = gsw.density.sigma1(SA, CT)
    
    return SA, CT, sigma0, T_f





#OLD FUNCTIONS
# def calc_friction_v(u_ice = np.array([]), v_ice = np.array([]),
#                     u_ocn = np.array([]), v_ocn = np.array([]),
#                     rho_o = 1023 * units('kg/m3'), Cio = 0.0055):

#     # ice-ocean velocity difference
#     U0_x = (u_ice - u_ocn)
#     U0_y = (v_ice - v_ocn)
    
#     # magnitude of ice-ocean velocity difference
#     U0 = np.sqrt(U0_x**2 + U0_y**2)
    
# #     # I think the following only matters for direction of stress
# #     # they say (+) turning angle, but in oceanography (+) means clockwise - yes?
# #     # so to get CW turning, apply (-) algebraic rotation
# #     beta = -23 * np.pi/180 
    
# #     # angle of ice-ocean velocity difference
# #     theta = np.arctan2(U0_y, U0_x)
    
# #     # rotated ice-ocean velocity difference
# #     rot_x = U0 * np.cos(theta+beta)
# #     rot_y = U0 * np.sin(theta+beta)
    
# #     # magnitude of rotated ice-ocean velocity difference
# #     R0U0 = np.sqrt(rot_x**2 + rot_y**2)
    
#     # simplified: IGNORING DIRECTION:
#     R0U0 = U0
    
#     # Zhong description of parameterized friction velocity from ice-ocean drag:
#     # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021GL096216
#     us = np.sqrt(Cio * U0 * R0U0)
    
#     return us
    



