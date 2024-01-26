#////////////////////////////
#  get_MODISgeo_Corners  ///
#//////////////////////////
#---------------------------------------------------------------------
# Load corners of lat, lon arrays from MODIS geo (hdf) files.
#---------------------------------------------------------------------
#///////////////////////
#  colorscale_range ///
#/////////////////////
#---------------------------------------------------------------------
# Given list of loaded MODIS images, find min/max of reflectance/radiance
#---------------------------------------------------------------------
#////////////////////
#  get_hdf_data  ///
#//////////////////
#---------------------------------------------------------------------
# Load data from hdf file (can load single attribute or full dataset). 
#---------------------------------------------------------------------
#///////////////////////
#  load_MODISband   ///
#/////////////////////
#---------------------------------------------------------------------
# Load a band from MODIS imagery.
#---------------------------------------------------------------------
#/////////////////////////
#  crop_makepoly_geo  ///
#///////////////////////
#---------------------------------------------------------------------
# Read in geofile lats and lons from MODIS imagery.
#---------------------------------------------------------------------
#////////////////////
#  get_MODISgeo  ///
#//////////////////
#---------------------------------------------------------------------
# Load lat, lon arrays from MODIS geo (hdf) files.
#---------------------------------------------------------------------
#/////////////////////
#  get_MODISdate  ///
#///////////////////
#---------------------------------------------------------------------
# Grab date from MODIS filename, create datetime object.
#---------------------------------------------------------------------


#////////////////////////////
#  get_MODISgeo_Corners  ///
#//////////////////////////
#---------------------------------------------------------------------
# Load corners of lat, lon arrays from MODIS geo (hdf) files.
#---------------------------------------------------------------------
# DEPENDENCIES
from pyhdf.SD import SD, SDC
#---------------------------------------------------------------------
def get_MODISgeo_Corners(geofile):

    """Load corners of lat, lon arrays from MODIS geo (hdf) files.
    Reads in Terra/MODIS (MOD03) or Aqua/MODIS (MOD03) geo files, 
    may also work for hdf geolocation files from other satellites
    returns longitudes in range (0, 360).
    
INPUT:
- geofile: filename with directory 
            (e.g. '/Users/kenzie/MOD03.A2000059.1745.061.2017171195808.hdf')

OUTPUT:
- corner1: [lat1, lon1]
- corner2: [lat2, lon2]
- corner3: [lat3, lon3]
- corner4: [lat4, lon4]

DEPENDENCIES:
from pyhdf.SD import SD, SDC

Latest recorded update:
06-03-2022

    """
    # open geo file
    #--------------   
    try:
        f = SD(geofile,SDC.READ) 
    # raise an error if file can't be opened 
    except Exception as e:
        print(e, ", error opening ", geofile, sep='')

    # grab each corner  [lat, lon]
    #-----------------------------
    corner1 = [f.select(0)[0,0],f.select(1)[0,0]]
    corner2 = [f.select(0)[0,-1],f.select(1)[0,-1]]
    corner3 = [f.select(0)[-1,-1],f.select(1)[-1,-1]]
    corner4 = [f.select(0)[-1,0],f.select(1)[-1,0]]

    # make lon range 0 < lon < 360
    #-----------------------------
    if corner1[1] < 0:
        corner1[1] += 360
    if corner2[1] < 0:
        corner2[1] += 360
    if corner3[1] < 0:
        corner3[1] += 360
    if corner4[1] < 0:
        corner4[1] += 360
    
    # close geo file and return corners
    #----------------------------------
    f.end() 
    return corner1, corner2, corner3, corner4

#///////////////////////
#  colorscale_range ///
#/////////////////////
#---------------------------------------------------------------------
# Given list of loaded MODIS images, find min/max of reflectance/radiance
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
#---------------------------------------------------------------------
def colorscale_range(image_list):
    
    """Given list of loaded MODIS images, find minimum/maximum 
    of reflectance/radiance values from all images in list.
    (this can be used to normalize image color scale on plot 
    containing multiple overlain images)
    
INPUT:
- image_list: list of loaded MODIS imagery files=
              (each file contains reflectance/radiance values to be plotted)

OUTPUT:
- _min, _max: minimum, maximum ref or rad values of all images in list

DEPENDENCIES:
import numpy as np

Latest recorded update:
06-03-2022

    """
    
    # initialize dummy values for min and max radiance values
    _min, _max = np.inf, 0
    # for each file in list
    for ii in range(0,len(image_list)):
        # find the minimum, maximum radiance values
        minimum, maximum = np.min(image_list[ii]), np.max(image_list[ii])
        # if radiance values are outside range of previous image file min/max, 
        # extend min/max to the wider range
        if minimum < _min:
            _min = minimum
        if maximum > _max:
            _max = maximum
                
    return _min, _max


#////////////////////
#  get_hdf_data  ///
#//////////////////
#---------------------------------------------------------------------
# Load data from hdf file (can load single attribute or full dataset). 
#---------------------------------------------------------------------
# DEPENDENCIES
from pyhdf.SD import SD, SDC
#---------------------------------------------------------------------
def get_hdf_data(file,dataset,attr):
    
    """Load data from hdf file (can load single attribute or full dataset).
    
INPUT:
- file: hdf filename with directory  
        (e.g. '/Users/kenzie/MOD021KM.A2000066.2255.061.2017171220013.hdf')
- dataset: desired data set within HDF file 
           (e.g. 'EV_250_Aggr1km_RefSB')
- attr: None OR desired attribute within dataset 
        (e.g. None, 'reflectance_scales')

OUTPUT:
- specified dataset or attribute

DEPENDENCIES:
from pyhdf.SD import SD, SDC

Latest recorded update:
06-03-2022

    """
    
    f = SD(file,SDC.READ)
    data = f.select(dataset)
    # if no attribute, grab full data
    if attr == None:            
        data_or_attr = data[:]
    # or grab attribute from full data
    else:                       
        index = data.attr(attr).index()
        data_or_attr = data.attr(index).get()
    f.end()
    
    return data_or_attr



#///////////////////////
#  load_MODISband   ///
#/////////////////////
#---------------------------------------------------------------------
# Load a band from MODIS imagery.
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
# homemade: get_hdf_data
#---------------------------------------------------------------------

def load_MODISband(file, dataset, band, refrad):

    """Load a band from MODIS imagery. Applies scale factor and offsets,
    makes mask for invalid/missing data values.
    
INPUT:
- file: filename with directory 
        (e.g. '/Users/kenzie/MOD021KM.A2000066.2255.061.2017171220013.hdf')
- dataset: desired data set within HDF file 
           (e.g. 'EV_250_Aggr1km_RefSB')
- band: band number formatted as string 
        (e.g. '30')
- refrad: reflectance or radiance datatype
          ('reflectance' or 'radiance')

OUTPUT:
- band_data: ref or rad band data, "bad" data masked

DEPENDENCIES:
import numpy as np
# homemade: get_hdf_data

Latest recorded update:
06-03-2022

    """
    
    # import data
    #---------------------------------------------------------------------
    band_names = get_hdf_data(file, dataset, 'band_names')
    band_data = get_hdf_data(file, dataset, None)[band_names.split(",").index(band), :, :].astype(np.double)
    if refrad == 'reflectance':
        scales = get_hdf_data(file, dataset, 'reflectance_scales')[band_names.split(",").index(band)]
        offsets = get_hdf_data(file, dataset, 'reflectance_offsets')[band_names.split(",").index(band)]
    elif refrad == 'radiance':
        scales = get_hdf_data(file, dataset, 'radiance_scales')[band_names.split(",").index(band)]
        offsets = get_hdf_data(file, dataset, 'radiance_offsets')[band_names.split(",").index(band)]   
    else:
        print('REFLECTANCE OR RADIANCE NOT SPECIFIED')
        
    validmin = get_hdf_data(file, dataset, 'valid_range')[0]
    validmax = get_hdf_data(file, dataset, 'valid_range')[1]
    fillval = get_hdf_data(file, dataset, '_FillValue')
    
    # make mask of data, elminating "bad" data  
    # ---------------------------------------------------------------------
    # identify fill values or values outside valid range
    invalid = np.logical_or(band_data > validmax, band_data < validmin)
    invalid = np.logical_or(invalid, band_data == fillval)
    # replace these ^ with NaNs
    band_data[invalid] = np.nan
    # apply offset and scale
    band_data = (band_data - offsets) * scales 
    # make data a masked array to ignore NaNs in calculations
    band_data = np.ma.masked_array(band_data, np.isnan(band_data))
    
    return band_data


#/////////////////////////
#  crop_makepoly_geo  ///
#///////////////////////
#---------------------------------------------------------------------
# Read in geofile lats and lons from MODIS imagery.
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
from shapely.geometry.polygon import Polygon
#---------------------------------------------------------------------
def crop_makepoly_geo(geolat, geolon, alSt, alEn, crSt, crEn, make_polygon):
    
    """Read in geofile lats and lons from MODIS imagery. 
    Crop image and return lat/lon grid with polygon perimeter.
    Option to not return polygon perimeter if not desired.
    
INPUT:
- geolat: array of lat values
- geolon: array of lon values 
- alSt, alEn: first and last indices of along-track dimension of level1b image
- crSt, crEn: first and last indices of cross-track dimension of level1b image
- make_polygon: 'makepoly' to make polygon, or set to anything if not wanted

OUTPUT:
- LAT, LON: Lat/lon grids from MODIS geolocation file, cropped to image dimensions
- poly: list of 8 coordinates along image perimeter
- polygon: polygon geoshape of image perimeter

DEPENDENCIES:
import numpy as np
from shapely.geometry.polygon import Polygon

Latest recorded update:
06-03-2022

    """
    
    # crop image based off coords
    #-----------------------------
    LAT = geolat[alSt:alEn, crSt:crEn]
    LON = geolon[alSt:alEn, crSt:crEn]
    
    # make polygon
    #-----------------------------
    if make_polygon == 'makepoly':
        # make lat/lon list of 8 points (corners, sides) along perimeter of image
        poly = np.array([[LON[-1, -1], LAT[-1, -1]], 
                                   [LON[round(2*LAT.shape[0]/3), -1], LAT[round(2*LAT.shape[0]/3), -1]],
                                   [LON[round(LAT.shape[0]/3), -1], LAT[round(LAT.shape[0]/3), -1]],
                                   [LON[0, -1], LAT[0, -1]],
                                   [LON[0, round(2*LAT.shape[1]/3)], LAT[0, round(2*LAT.shape[1]/3)]],
                                   [LON[0, round(LAT.shape[1]/3)], LAT[0, round(LAT.shape[1]/3)]],
                                   [LON[0, 0], LAT[0, 0]], 
                                   [LON[round(LAT.shape[0]/3), 0], LAT[round(LAT.shape[0]/3), 0]],
                                   [LON[2*round(LAT.shape[0]/3), 0], LAT[2*round(LAT.shape[0]/3), 0]],
                                   [LON[-1, 0], LAT[-1, 0]],
                                   [LON[-1, round(LAT.shape[1]/3)], LAT[-1, round(LAT.shape[1]/3)]],
                                   [LON[-1, round(2*LAT.shape[1]/3)], LAT[-1, round(2*LAT.shape[1]/3)]]
                                  ])
        # make polygon geoshape with specified poly coordinates 
        polygon = Polygon(poly)
        
        return LAT, LON, poly, polygon
    
    else:
        
        return LAT, LON
    
    
#////////////////////
#  get_MODISgeo  ///
#//////////////////
#---------------------------------------------------------------------
# Load lat, lon arrays from MODIS geo (hdf) files.
#---------------------------------------------------------------------
# DEPENDENCIES
from pyhdf.SD import SD, SDC
#---------------------------------------------------------------------
def get_MODISgeo(geofile):

    """Load lat, lon arrays from MODIS geo (hdf) files.
    Reads in Terra/MODIS (MOD03) or Aqua/MODIS (MOD03) geo files,
    May also work for hdf geolocation files from other satellites.
    Returns longitudes in range (0, 360).
    
INPUT:
- geofile: filename with directory  
           (e.g. '/Users/kenzie/MOD03.A2000059.1745.061.2017171195808.hdf')

OUTPUT:
- geolat: array of lat values
- geolon: array of lon values

DEPENDENCIES:
from pyhdf.SD import SD, SDC

Latest recorded update:
06-03-2022

    """
    
    # open geo file
    #--------------   
    try:
        f = SD(geofile,SDC.READ) 
    # raise an error if file can't be opened 
    except Exception as e:
        print(e, ", error opening ", geofile, sep='')

    # open geo file
    #--------------
    geolat = f.select(0)[:]       # pull out lat    
    geolon = f.select(1)[:]       # pull out lon 
    
    # make all longitudes range (0,360)
    #----------------------------------
    for i in range(0,geolon.shape[0]): 
         for j in range(0,geolon.shape[1]):
                if geolon[i,j] < 0:
                    geolon[i,j] = 360 + geolon[i,j]  
                    
    # close geo file and return lat, lon
    #-----------------------------------
    f.end() 
    
    return geolat, geolon


#/////////////////////
#  get_MODISdate  ///
#///////////////////
#---------------------------------------------------------------------
# Grab date from MODIS filename, create datetime object.
#---------------------------------------------------------------------
# DEPENDENCIES
import datetime as dt
from datetime import datetime
#---------------------------------------------------------------------
def get_MODISdate(filename):
    
    """Grab date from MODIS filename, create datetime object.
    MODIS filename can be from either geolocation or imagery files
    from level1b modis products, as long as they include date after
    '.A' in the filename.
    
INPUT:
- file: MODIS filename (without path) 
        (e.g. 'MOD021KM.A2006090.2150.061.2017263004124.hdf')

OUTPUT:
- imagedate: datetime object of MODIS image

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




# #/////////////////////////////
# #  plot_singleband_MODIS  ///
# #///////////////////////////
# #---------------------------------------------------------------------
# # Make plot of projected single-band MODIS image from imagery and geo file.
# #---------------------------------------------------------------------
# # DEPENDENCIES:
# #-------------
# import numpy as np
# import cartopy
# import cartopy.crs as ccrs
# import matplotlib as mpl
# from matplotlib import pyplot as plt
# import matplotlib.colors                                
# #---------------------------------------------------------------------

# def plot_singleband_MODIS(LAT, LON, _image_, 
#               TargetCCRS = ccrs.NorthPolarStereo(central_longitude=215), 
#               SourceCCRS = ccrs.PlateCarree(), 
#               ImageSize = [10,12], Extent = None, 
#               cmap = 'Greys',  cscale = [1.8, 6.8], shading = 'gouraud',zorder=1): 
    
    
#     """Make plot of projected single-band MODIS image from imagery and geo file.
    
# INPUT:
# - LAT, LON: lat/lon grids of image to plot (or list of grids if overlaying multiple images simultaneously)
# - _image_ : image to plot (or list of images if overlaying multiple images simultaneously)
# - TargetCCRS: Desired CRS of plot (https://scitools.org.uk/cartopy/docs/latest/crs/projections.html)
#                 (default: ccrs.NorthPolarStereo(central_longitude=215))
                
# - SourceCCRS: CRS of input data, should use ccrs.PlateCarree() (default: ccrs.PlateCarree())
# - ImageSize: figure size of image (default: [10, 12])
# - Extent: plot extent ([lonmin, lonmax, latmin, latmax]) or None if map should
#           not be cropped to specific extent (default: None)
#           (units of SourceCCRS, here lat/lon like ccrs.PlateCarree())
# - cmap: colormap t0 use (default: 'Greys')
# - cscale: normalization for data in colormap (default: [min, max] = [1.8, 6.8])
#             Mx1 lists of floats (M=2,3,4) either as:
#               - [vmin, vmax] 
#               - [vmin, midpoint, vmax] 
#               - [vmin, midpoint1, midpoint2, vmax]

# - shading: pcolormesh shading style (default: 'gouraud')
# - zorder: drawing order of layer (default: 1)


# OUTPUT:
# - fig, ax: figure and axis with projected single-band MODIS image

# DEPENDENCIES:
# import numpy as np
# import cartopy
# import cartopy.crs as ccrs
# import matplotlib as mpl
# from matplotlib import pyplot as plt
# import matplotlib.colors

# Latest recorded update:
# 06-03-2022

#     """
    
#     # create figure
#     #--------------
# #     fig, ax = plt.subplots(subplot_kw=dict(projection=TargetCCRS), figsize=(ImageSize[0], ImageSize[1]))

#     # set extent
#     #-----------
# #     if Extent != None:
# #         ax.set_extent(Extent, crs=SourceCCRS)

#     # colormap normalization
#     #-----------------------
#     # class for colormap normalizing scaling options
#     class TwopointNormalize(matplotlib.colors.Normalize):
#         def __init__(self, vmin=None, vmax=None, vmid1=None, vmid2=None, clip=False):
#             self.vmid1 = vmid1
#             self.vmid2 = vmid2
#             super().__init__(vmin, vmax, clip)
#         def __call__(self, value, clip=None):
#             x, y = [self.vmin, self.vmid1, self.vmid2, self.vmax], [0, 0.33,0.66, 1]
#             return np.ma.masked_array(np.interp(value, x, y))

#     # determine color normalization
#     if len(cscale)==2:
#         norm=matplotlib.colors.Normalize(vmin=cscale[0], vmax=cscale[1])
#     elif len(cscale)==3:
#         norm=matplotlib.colors.TwoSlopeNorm(vmin=cscale[0], vcenter=cscale[1],vmax=cscale[2])
#     else:
#         norm = TwopointNormalize(vmin=cscale[0], 
#                                         vmid1=cscale[1], 
#                                         vmid2=cscale[2], 
#                                         vmax=cscale[3])
        
#     # plot level1b image
#     #-------------------
    
#     # for each file in _image_ list
#     for ii in range(0,len(_image_)):
#         ax.pcolormesh(LON[ii], LAT[ii], _image_[ii], 
#                       cmap = cmap, norm = norm,
#                       transform=SourceCCRS, zorder = zorder, shading = shading)
#         # use 'gouraud' shading in pcolormesh since it doesn't like that it's given cell centers rather than corners
        
#     # remove automatic image border
#     ax.outline_patch.set_linewidth(0)
    
# #     return fig, ax


