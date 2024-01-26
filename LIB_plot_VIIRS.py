#////////////////////////////
#  get_VIIRSgeo_Corners  ///
#//////////////////////////
#---------------------------------------------------------------------
# Load corners of lat, lon arrays from VIIRS geo (nc) files.
#---------------------------------------------------------------------
#///////////////////////////
#  load_VIIRS_band   ///
#/////////////////////////
#---------------------------------------------------------------------
# Load a band from VIIRS imagery.
#---------------------------------------------------------------------
#////////////////////////
#  get_VIIRS_geo  ///
#//////////////////////
#---------------------------------------------------------------------
# Load lat, lon arrays from VIIRS geo (nc) files.
#---------------------------------------------------------------------
#/////////////////////////
#  get_VIIRS_date  ///
#///////////////////////
#---------------------------------------------------------------------
# Grab date from VIIRS filename, create datetime object.
#---------------------------------------------------------------------




#////////////////////////////
#  get_VIIRSgeo_Corners  ///
#//////////////////////////
#---------------------------------------------------------------------
# Load corners of lat, lon arrays from VIIRS geo (nc) files.
#---------------------------------------------------------------------
# DEPENDENCIES
import netCDF4
#---------------------------------------------------------------------
def get_VIIRSgeo_Corners(geofile):

    """Load corners of lat, lon arrays from VIIRS geo (nc) files.
    Returns longitudes in range (0, 360).
    
INPUT:
- geofile: filename with directory 
            (e.g. '/Users/kenzie/VNP03MOD.A2021078.1218.002.2021127035344.nc')

OUTPUT:
- corner1: [lat1, lon1]
- corner2: [lat2, lon2]
- corner3: [lat3, lon3]
- corner4: [lat4, lon4]

DEPENDENCIES:
import netCDF4

Latest recorded update:
09-29-2022

    """
    
    # open geo file
    #--------------   
    try:
        # open geo file
        f = netCDF4.Dataset(geofile)
    # raise an error if file can't be opened 
    except Exception as e:
        print(e, ", error opening ", geofile, sep='')

    # open geo file
    #--------------
    lats = f.groups['geolocation_data'].variables['latitude'][:]
    lons = f.groups['geolocation_data'].variables['longitude'][:]
    lons[lons<0]+=360

    # grab each corner  [lat, lon]
    #-----------------------------
    corner1 = [ lats[0,0],   lons[0,0]   ]
    corner2 = [ lats[0,-1],  lons[0,-1]  ]
    corner3 = [ lats[-1,-1], lons[-1,-1] ]
    corner4 = [ lats[-1,0],  lons[-1,0]  ]
    
    # close geo file and return corners
    #----------------------------------
    f.close() 
    
    return corner1, corner2, corner3, corner4





#////////////////////////
#  load_VIIRS_band   ///
#//////////////////////
#---------------------------------------------------------------------
# Load a band from VIIRS imagery.
#---------------------------------------------------------------------
# DEPENDENCIES
import netCDF4
#---------------------------------------------------------------------

def load_VIIRS_band(file, band = 'M15'):

    """Load a band from VIIRS imagery file. Appears that scale factor and offsets, 
    and mask for invalid/missing data values are automatically applied. 
    
INPUT:
- file: filename with directory 
        (e.g. '/Users/kenzie/VNP02MOD.A2021078.1218.002.2021128180808.nc')
- band: band number formatted as string (default: 'M15')
        (e.g. 'M15')

OUTPUT:
- VIIRSimg: radiances data of specified, "bad" data masked

DEPENDENCIES:
import netCDF4

Latest recorded update:
09-29-2022

    """
    # open geo file
    #--------------   
    try:
        # open geo file
        f = netCDF4.Dataset(file)
    # raise an error if file can't be opened 
    except Exception as e:
        print(e, ", error opening ", file, sep='')
        
    # appears data is already masked and scaled...
    VIIRSimg = f.groups['observation_data'].variables[band][:]
    f.close()
 
    return VIIRSimg



#/////////////////////
#  get_VIIRS_geo  ///
#///////////////////
#---------------------------------------------------------------------
# Load lat, lon arrays from VIIRS geo (nc) files.
#---------------------------------------------------------------------
# DEPENDENCIES
import netCDF4
#---------------------------------------------------------------------
def get_VIIRS_geo(geofile):

    """Load lat, lon arrays from VIIRS geo (nc) files.
    Returns longitudes in range (0, 360).
    
INPUT:
- geofile: filename with directory  
           (e.g. '/Users/kenzie/VNP03MOD.A2021078.1218.002.2021127035344.nc')

OUTPUT:
- lats: array of lat values
- lons: array of lon values

DEPENDENCIES:
import netCDF4

Latest recorded update:
09-29-2022

    """
    
    # open geo file
    #--------------   
    try:
        # open geo file
        f = netCDF4.Dataset(geofile)
    # raise an error if file can't be opened 
    except Exception as e:
        print(e, ", error opening ", geofile, sep='')

    # open geo file
    #--------------
    lats = f.groups['geolocation_data'].variables['latitude'][:]
    lons = f.groups['geolocation_data'].variables['longitude'][:]
    lons[lons<0]+=360

    # close geo file and return lat, lon
    #-----------------------------------
    f.close()
    
    return lats, lons





#//////////////////////
#  get_VIIRS_date  ///
#////////////////////
#---------------------------------------------------------------------
# Grab date from VIIRS filename, create datetime object.
#---------------------------------------------------------------------
# DEPENDENCIES
import datetime as dt
from datetime import datetime
#---------------------------------------------------------------------
def get_VIIRS_date(filename):
    
    """Grab date from VIIRS filename, create datetime object.
    VIIRS filename can be from either geolocation or imagery files, 
    as long as they include date after
    '.A' in the filename.
    
INPUT:
- file: VIIRS filename (without path) 
        (e.g. 'VNP02MOD.A2021078.1218.002.2021128180808.nc')

OUTPUT:
- imagedate: datetime object of VIIRS image


DEPENDENCIES:
import datetime as dt
from datetime import datetime

Latest recorded update:
09-29-2022

    """
    
    # find beginning of date in filename
    #------------------------------------
    di = filename.index('.A')+2
    # grab year, day, hour, and minutes of image acquisition
    #---------------------------------------------------------
    YYYY = str(filename[di:di+4])
    DDD = str(filename[di+4:di+7])
    HH = str(filename[di+8:di+10])
    MM = str(filename[di+10:di+12])
    # create and return datetime object
    #----------------------------------
    imagedate = dt.datetime.strptime(YYYY+' '+DDD+' '+HH+' '+MM, '%Y %j %H %M')
    
    return imagedate