# DEPENDENCIES:
import os
import numpy as np
import numpy.ma as ma
import datetime 
import xarray as xr
import pandas as pd

def grab_iceage(date = datetime.datetime(year = 2000, month = 1, day = 1),
                select_date = 'nearest',
                iceage_datapath = '/Volumes/Jewell_EasyStore/NSIDC-0611_seaice_age/', 
                iceage_filename = 'iceage_nh_12.5km_{}0101_{}1231_v4.1.nc',
                return_vars = ['lon', 'lat', 'age', 'ds', 'age-key', 'selected_date'], 
                mask_flags = True,
                quiet = True):
    """Grab sea ice age data from locally-stored NSIDC data (NSIDC-0611, doi: 10.5067/UTAV7490FEPB)
        COULD BE FIXED LATER: DON'T JUST OPEN DATA OF date's YEAR. THIS WILL NOT EFFECTIVELY GRAB NEAREST DATES AROUND START/END OF YEARS.

INPUT: 
- date: data of data to open
- select_date: way to select nearest time (of weekly data) to provided date:
    - 'nearest': grab nearest date to that provided (default)
    - 'before': grab nearest date before provided
    - 'after': grab nearest date after provided
    * if date is exact match to date in file, each option will return same date
- iceage_datapath: path to local data files
- iceage_filename: naming convention of stored data files, with {} indicating year as YYYY
- return_vars: variables/attributes to return in specified order. (list)
    Can include any or all OUTPUT variables in any order: 
    default: ['lon', 'lat', 'age', 'ds', 'age-key', 'selected_date']
- mask_flags: whether or not to mask flagged (age=20, land and age=21, near-land/lakes) (default: True)
- quiet: bool, whether or not to hide print statements (default: True)
    

OUTPUT:
List of any or all of variables specified in return_vars:
- lon: M x N array of longitudes (0 - 360)
- lat: M x N array of latitudes
- age: M x N array of sea ice ages
- age-key: string listing sea ice age key
- selected_date: actual date selected associated with ice age data 
- ds: xarray data frame containing opened data

DEPENDENCIES:
import os
import numpy as np
import numpy.ma as ma
import datetime 
import xarray as xr
import pandas as pd

Latest recorded update:
06-27-2023
    """
    
    # 0 Open water or < 15% sea ice concentration
    # Sea ice age; higher age estimates are not precise, so older ice, 5th-year (4-5 years
    # old) and above, are generally considered together
    # 1 = ice that is 0-1 years old (first-year ice)
    # 2 = ice that is 1-2 years old (second-year ice)
    # 3 = ice that is 2-3 years old (third-year ice)
    # …
    # 16 = ice that is 15-16 years old (16th-year ice)
    # 20 Designates the grid cell contains only land
    # 21 Designates grid cells that contain ocean for which ice age was not calculated
    
    # open file from year of event
    # FIX LATER !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ds = xr.open_dataset(iceage_datapath+iceage_filename.format(date.year, date.year))
    ds.close()
    
    # change time index form cftime to datetime
    datetimeindex = ds.indexes['time'].to_datetimeindex()
    datetimeindex
    ds['time'] = datetimeindex
    
    # finding nearest date in file
    #-------------------------------------------------
    # find time difference between dates
    times = pd.to_datetime(ds.time.values)
    time_diffs = times-date

    # grab nearest date
    if str(select_date) == 'nearest':
        date_index = np.nanargmin(np.abs(time_diffs))
        string = 'to'

    # grab nearest date before, raising error it not possible
    elif str(select_date) == 'before':
        if (times[0] - date).days > 0:
            raise ValueError(f'\nERROR: Cannot grab date before {date.date()} within time range [{times[0].date()}, {times[-1].date()}]')
        else:
            time_diffsecs = np.array(time_diffs.total_seconds())
            time_diffsecs[time_diffsecs > 0] = np.nan
            date_index = np.nanargmax(time_diffsecs)
        string = 'before'

    # grab nearest date after, raising error it not possible
    elif str(select_date) == 'after':
        if (date - times[-1]).days > 0:
            raise ValueError(f'\nERROR: Cannot grab date after {date.date()} within time range [{times[0].date()}, {times[-1].date()}]')
        else:
            time_diffsecs = np.array(time_diffs.total_seconds())
            time_diffsecs[time_diffsecs < 0] = np.nan
            date_index = np.nanargmin(time_diffsecs)
            string = 'after'
            
    # grab nearest date and index of nearest date
    nearest_date = times[date_index]
    if not quiet:
        print(f' - select nearest date {string} {date.date()} >>> {nearest_date.date()}')
    # warn if selected date more than a week form input date
    if np.abs((date - nearest_date).days) > 7:
        print(' >>> warning: selected date more than 7 days from input date')
    
    # grabbing data
    #-------------------------------------------------
    # selecgt date
    ds_spec = ds.sel(time=ds.time[date_index])

    # extract age and coordinates from data set
    seaice_age = ds_spec.age_of_sea_ice.values
    ice_lon = ds_spec.longitude.values
    ice_lon[ice_lon<0]+=360
    ice_lat = ds_spec.latitude.values
    
    # apply land/lake mask if desired
    if mask_flags:
        flags = ds.age_of_sea_ice.flag_values
        seaice_age = ma.masked_where(np.isin(seaice_age, flags),seaice_age)
        
    # store all output in dictionary
    vars_dict = {}
    vars_dict['lon'] = ice_lon
    vars_dict['lat'] = ice_lat
    vars_dict['age'] = seaice_age
    vars_dict['age-key'] = ds.age_of_sea_ice.comment + ". 0: Open water or < 15% sea ice concentration. Higher age estimates are not precise, so older ice, 5th-year (4-5 years old) and above, are generally considered together"
    vars_dict['ds'] = ds
    vars_dict['selected_date'] = nearest_date
    
    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    if len(return_data) == 1:
        return_data = return_data[0]
    return return_data





