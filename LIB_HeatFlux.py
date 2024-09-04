

from metpy.units import units
import numpy as np
from metpy.calc import relative_humidity_from_dewpoint, mixing_ratio_from_relative_humidity, vapor_pressure


# calc_ssh
# calc_lwd
# calc_lwu
# seasonal_Cio
# calc_vars_from_T_S

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
