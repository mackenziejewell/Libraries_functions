# DEPENDENCIES:
import os
import h5py
from datetime import datetime

def generate_TB_filename(date, quiet = False):
    
    """Generate name of locally-stored / re-written TB files from date. 

INPUT: 
- date: date of file to open
- quiet: bool, whether or not to suppress prints (default: False)

OUTPUT:
- data: dictionary storing TB data 

DEPENDENCIES:
from datetime import datetime

Latest recorded update:
05-21-2024
    """
    
    if date < datetime(2012, 10, 8):
        path = '/Volumes/Seagate_Jewell/KenzieStuff/AMSR_TB/AE_SI12/'
    else:
        path = '/Volumes/Seagate_Jewell/KenzieStuff/AMSR_TB/AU_SI12/'

    filename = f"AMSR_TB_L3_{date.strftime('%Y%m%d')}.h5"
    file = path + filename
    
    if not os.path.isfile(file):

        if not quiet:
            print(f'!!! File does not exist: {file}')
            
        file = False
            
    return file


def grab_TB_file_data(file, 
                      frequencies = ['18', '36', '89'], 
                      polarizations = ['V', 'H'],
                      overpasses = ['ASC', 'DSC', 'DAY'], 
                      include_SIC = True):

    """Grab from home-re-written files based on TB data from NSIDC 12.5 km product. Read into dictionary. 

INPUT: 
- file: name of locally-stored file
- frequencies: list of frequencies in GHz to read into dict (default: ['18', '36', '89'])
- polarizations: list of polarizations of TB data to read in (default: ['V', 'H'])
- overpasses: list of overpasses to store in dict (default: ['ASC', 'DSC', 'DAY'])
- include_SIC: whether or not to read in SIC data

OUTPUT:
- data: dictionary storing TB data 

DEPENDENCIES:
import os
import h5py

Latest recorded update:
05-21-2024
    """
    
    hf = h5py.File(file,'r')

    # dictionary to store data
    data = {}

    # file data, could be adjusted later but applies to locally-made files
    resolution = '12km'
    hemisphere = 'N'


    # read in TB data of desired freq, pol, overpass
    #-----------------------------------------------
    for f in frequencies:
        data[f] = {}

        for p in polarizations:
            data[f][p] = {}

            for o in overpasses:
                
                variable = f'SI_{resolution}_{hemisphere}H_{f}{p}_{o}'

                data[f][p][o] = hf[variable][:] * 0.1 # scale factor 0.1
                # according to documentation:
                # https://nsidc.org/sites/default/files/ae_si6-v003-userguide.pdf
                # _FillValue   = 0
                # scale_factor = 0.1
                # add_offset   = 0

    # read in SIC if desired
    #-----------------------
    if include_SIC:    
        
        variable = f'SI_{resolution}_{hemisphere}H_ICECON_DAY'
        data['SIC'] = hf[variable][:]
    
    return data


def gradient_ratio(data, f1 = '89', f2 = '18', p = 'V', overpass = 'DAY'):
    
    """Calculation TB gradient ratio: GR(f1P, f2P) = (f1P - f2P) / (f1P + f2P)

    INPUT: 
    - data: dictionary storing TB data (created using "grab_TB_file_data")
    - f1: string of 1st TB frequency (GHz) to use in GR (default: '89')
    - f2: string of 2nd TB frequency (GHz) to use in GR (default: '18')
    - p: polarization to use in GR (default: 'V')
    - overpass: string of overpass to calculate GR (default: 'DAY')

    OUTPUT:
    - data: calculated GR

    DEPENDENCIES:

    Latest recorded update:
    05-21-2024
    """
    
    GR = (data[f1][p][overpass] - data[f2][p][overpass]) / (data[f1][p][overpass] + data[f2][p][overpass])
    
    return GR




def polarization_ratio(data, f = '18', overpass = 'DAY'):

    """Calculation TB polarization ratio: PR(f) = (fV - fH) / (fV + fH)

    INPUT: 
    - data: dictionary storing TB data (created using "grab_TB_file_data")
    - f: string of TB frequency (GHz) to calculate PR (default: '18')
    - overpass: string of overpass to calculate PR (default: 'DAY')

    OUTPUT:
    - data: calculated PR

    DEPENDENCIES:

    Latest recorded update:
    05-21-2024
    """
     
    PR = (data[f]['V'][overpass] - data[f]['H'][overpass]) / (data[f]['V'][overpass] + data[f]['H'][overpass])

    return PR


