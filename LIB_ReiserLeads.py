import os
import numpy.ma as ma
import numpy as np
import xarray as xr

def grab_ReiserLeadFreq(file = '/Volumes/Jewell_EasyStore/Reiser_LeadFrequency_20022019/Arctic_Relleads_2002_2019.nc',
                        return_vars = ['lon', 'lat', 'freq', 'ds']):
    
    """Function to grab lead frequencies from locally stored Reiser et al. (2020) data. 
    
    Corresponds to data publication on PANGEA:
    Reiser, Fabian; Willmes, Sascha; Heinemann, Günther (2020): Daily sea ice lead data for Arctic and Antarctic. PANGAEA, https://doi.org/10.1594/PANGAEA.917588

INPUT: 
- file: path to Reiser et al. lead data
    deafult: '/Volumes/Jewell_EasyStore/Reiser_LeadFrequency_20022019/Arctic_Relleads_2002_2019.nc'
- return_vars: list of variables/attributes to return in specified order
    Can include any or all OUTPUT variables in any order: 
    default: ['lon', 'lat', 'freq', 'ds']
    - 'lon' : (M x N array of longitudes)
    - 'lat' : (M x N array of latitudes)
    - 'freq' : (M x N array of relative lead frequencies)
    - 'ds' : xarray dataset
    
OUTPUT:
- return_data: list of all data variables corresponding to names specified in return_vars
    for default return_vars specified above, read all output into list of variables for example as [lat, lon, freq, ds]

example:
[lat, lon, freq, ds] = grab_ReiserLeadFreq(return_vars = ['lon', 'lat', 'freq', 'ds'])

DEPENDENCIES:
import os
import numpy.ma as ma
import numpy as np
import xarray as xr

Latest recorded update:
06-20-2023
 """
    
    assert type(file) == str, f"file should be string, not {type(file)}"
    assert os.path.isfile(file), f"file = '{file}' not found"
    
    
    # open file
    ds = xr.open_dataset(file)
    ds.close()
    
    # store all output in dictionary
    vars_dict = {}
    
    # though data says scale factor = 1, the data are actually stored 0 - 255, 
    # so must be divided by 255 to get lead frequency as fractions between 0 and 1.
    # 255 corresponds to land, so mask these values out
    vars_dict['freq'] =  ma.masked_where(ds.LeadFrequency.values==255,ds.LeadFrequency.values/255)
    vars_dict['lat'] = ds.Latitude.values
    vars_dict['lon'] = ds.Longitude.values
    vars_dict['ds'] = ds
    
    # save specified variables to list for output
    return_data = [vars_dict[var] for var in return_vars]
    if len(return_data) == 1:
        return_data = return_data[0]
    
    return return_data
           