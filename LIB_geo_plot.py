#////////////////////////////
#  fix_cartopy_vectors   ///
#//////////////////////////
#---------------------------------------------------------------------
# function to output vector components for plotting in cartopy. 
#---------------------------------------------------------------------
#////////////////
#  add_land  ///
#//////////////
#---------------------------------------------------------------------
# Add land feature to cartopy figure
#---------------------------------------------------------------------
#/////////////////
#  add_coast  ///
#///////////////
#---------------------------------------------------------------------
# Add coast feature to cartopy figure
#---------------------------------------------------------------------
#////////////////
#  add_grid  ///
#//////////////
#---------------------------------------------------------------------
# Add specified gridlines to cartopy figure
#---------------------------------------------------------------------
#////////////////
#  add_date  ///
#//////////////
#---------------------------------------------------------------------
# Add date label to cartopy plot.
#---------------------------------------------------------------------
#///////////////
#  add_mslp ///
#/////////////
#---------------------------------------------------------------------
# Add ECMWF mslp contours to cartopy plot.
#---------------------------------------------------------------------
#///////////////
#  add_vectors ///
#/////////////
#---------------------------------------------------------------------
# Add vectors to cartopy plot.
#---------------------------------------------------------------------
#////////////////////////
#  plot_scalar_mesh  ///
#////////////////////// 
#---------------------------------------------------------------------
# function for plotting regular or geographic data with pcolormesh
#---------------------------------------------------------------------
#/////////////////////
#  plot_contours  ///
#///////////////////
#---------------------------------------------------------------------
# Function for plotting regular or geographic data contours
#---------------------------------------------------------------------
#//////////////////
#  add_LFice   ///
#////////////////
#---------------------------------------------------------------------
# Plot 'landfast ice' perimeter around Alaskan coast from bathymetry data
#---------------------------------------------------------------------


#////////////////////////////
#  fix_cartopy_vectors   ///
#//////////////////////////
#---------------------------------------------------------------------
# function to output vector components for plotting in cartopy. 
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import numpy.ma as ma
#---------------------------------------------------------------------
def fix_cartopy_vectors(u, v, uvlats):
    
    """Function to output vector components for plotting in cartopy. 
    
    Reads in vectors and associated latitudes, return fixed vectors. 
    
    Cartopy doesn't know meters per degree increase toward the pole 
    in zonal direction as cosine(lat) when reprojecting for vectors 
    given as m/s (where  a meter covers different amount of a degree 
    depending on latitude), we need to rescale the u (east-west) speeds. 
    otherwise at high latitudes, u will be drawn much larger than 
    they should which will give incorrect angles

INPUT: 
- u: (N x M) array of eastward velocity component (m/s)
- v: (N x M) array of northward velocity component (m/s)
- uvlats: (N x M) array latitudes associated with u,v vectors

OUTPUT:
- u_fixed: (N x M) array of u with correct angle and magnitude for plotting
- v_fixed: (N x M) array of v with correct angle and magnitude for plotting


DEPENDENCIES:
import numpy as np
import numpy.ma as ma

Latest recorded update:
04-20-2022
    """
    
    # FIX ANGLE
    # fix u scale to be correct relative to v scale
    #----------------------------------------------
    # for u (m/s), m/degree decreases as cos(lat)
    u_fixed = u/np.cos(uvlats/180*np.pi) 
    v_fixed = v  # v does not change with latitude

    # FIX MAGNITUDE
    # scale fixed u,v to have correct magnitude 
    #-----------------------
    # original magnitude 
    orig_mag = ma.sqrt(u**2+v**2)
    # new magnitude
    fixed_mag = ma.sqrt(u_fixed**2+v_fixed**2)
    u_fixed = u_fixed*(orig_mag/fixed_mag)
    v_fixed = v_fixed*(orig_mag/fixed_mag)

    return u_fixed, v_fixed


#////////////////
#  add_land  ///
#//////////////
#---------------------------------------------------------------------
# Add land features to cartopy figure
#---------------------------------------------------------------------
# DEPENDENCIES:
#-------------
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
from shapely import wkt
#---------------------------------------------------------------------
def add_land(ax, scale = '50m', color='gray', alpha=1, fill_dateline_gap = True, zorder=2):
    
    """Add land feature to cartopy figure
    
INPUT:
- ax: cartopy figure axis
- scale = NaturalEarthFeature land feature scale (e.g. '10m', '50m', '110m')
        (default: '50m')
- color = land color (e.g. 'k' or [0.9,0.6,0.5]) (default: 'gray')
- alpha = land opacity (default: 1)
- zorder: drawing order of land layer (default: 2)
- fill_dateline_gap: specify whether to fill gap in cartopy land feature along 
   dateline that crosses Russia and Wrangel Island (default: True)

OUTPUT:
- input plot with added land layer

DEPENDENCIES:
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib as mpl
from matplotlib import pyplot as plt
from shapely import wkt

Latest recorded update:
05-11-2022

    """
    
        
    # grab land from cfeat.NaturalEarthFeature
    #-----------------------------------------
    ax.add_feature(cfeat.NaturalEarthFeature(category='physical', name='land', 
                                             scale=scale, facecolor=color),
                                             alpha = alpha, zorder = zorder)

    # if specified, fill dateline gap in land feature with shapely polygons
    if fill_dateline_gap == True:
        # generate polygon to fill line across Wrangel Island and line across Russia
        WKT_fill_Wrangel = 'POLYGON ((-180.1 71.51,-180.1 71.01,-179.9 71.01,-179.9 71.51,-180.1 71.51))'
        poly1 = wkt.loads(WKT_fill_Wrangel)
        ax.add_geometries([poly1], crs=ccrs.PlateCarree(), 
              facecolor=color, edgecolor=color, alpha = alpha, zorder=zorder)
        WKT_fill_Russia = 'POLYGON ((-180.1 65.1,-180.1 68.96,-179.9 68.96,-179.9 65.1,-180.1 65.1))'
        poly2 = wkt.loads(WKT_fill_Russia)
        ax.add_geometries([poly2], crs=ccrs.PlateCarree(), 
              facecolor=color, edgecolor=color, alpha = alpha, zorder=zorder)


        
#/////////////////
#  add_coast  ///
#///////////////
#---------------------------------------------------------------------
# Add coast feature to cartopy figure
#---------------------------------------------------------------------
# DEPENDENCIES:
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
#---------------------------------------------------------------------
def add_coast(ax, scale = '50m', color='gray', linewidth = 1, alpha=1, zorder=3):

    """Add land feature to cartopy figure
    
INPUT:
- ax: cartopy figure axis
- scale = NaturalEarthFeature coast feature scale (e.g. '10m', '50m', '110m')
        (default: '50m')
- color = coastline color (e.g. 'k' or [0.9,0.6,0.5]) (default: 'gray')
- linewidth = coastline linewidth (default: 1)
- alpha = coastline opacity (default: 1)
- zorder: drawing order of coast layer (default: 3)

OUTPUT:
- input plot with added coast layer

DEPENDENCIES:
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib as mpl
from matplotlib import pyplot as plt
from shapely import wkt

Latest recorded update:
05-06-2022

    """

    # coastline
    #----------
    ax.coastlines(scale, color=color, linewidth=linewidth, alpha = alpha, zorder = zorder)



#////////////////
#  add_grid  ///
#//////////////
#---------------------------------------------------------------------
# Add specified gridlines to cartopy figure
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.ticker as mticker
import matplotlib as mpl
from matplotlib import pyplot as plt
#---------------------------------------------------------------------

def add_grid(ax, lats = None, lons = None, linewidth = 1, color = 'gray', alpha=0.5, zorder = 4): 
    
    """Add specified gridlines to cartopy figure.
    
INPUT:
- ax: cartopy figure axis
- lats: None or array of latitudes to plot lines (default: None)
- lons: None or array of latitudes to plot lines (default: None)
- linewdith: grid line linewidths (default: 1)
- color: grid line color (default: 'gray')
- alpha: line transparency (default: 0.5)
- zorder: drawing order of gridlines layer (default: 4)

OUTPUT:
- input plot with added grid

DEPENDENCIES:
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.ticker as mticker
import matplotlib as mpl
from matplotlib import pyplot as plt

Latest recorded update:
04-22-2022

    """
        
    # give gridline specifications
    #-----------------------------
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=linewidth, color=color, alpha=alpha, zorder = zorder)

    # add the longitude gridlines
    #----------------------------
    if lons is None:
        gl.xlocator = mticker.FixedLocator([])
    else:
        # shift all longitudes from [0,360] to [180,-180]
        lons = np.concatenate((lons[(lons>180)]-360,lons[(lons<=180)]))
        gl.xlocator = mticker.FixedLocator(lons)

        
    # add the latitude gridlines
    #----------------------------
    if lats is None:
        gl.ylocator = mticker.FixedLocator([])
    else:
        gl.ylocator = mticker.FixedLocator(lats)
        



        

#////////////////
#  add_date  ///
#//////////////
#---------------------------------------------------------------------
# Add date label to cartopy plot.
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
from datetime import datetime
from matplotlib.offsetbox import AnchoredText
#---------------------------------------------------------------------

def add_date(fig, ax, dt_obj, date_format = '%b %d, %Y (%H:%M UTC)', method = 'anchor',
             boxstyle="round,pad=0.,rounding_size=0.2", facecolor = 'black', edgecolor = 'black',
             zorder = 10,
             
             anchor_loc = 4, anchor_prop = {'size': 20, 'color':'white'},
             

             x = 0.02, y= 0.05, textcolor = 'white',fontsize=15): 
    
    """Add date label to cartopy plot.
    
INPUT:
- fig: cartopy figure
- ax: cartopy figure axis
- dt_obj: datetime object of date for plotted data 
            OR
          string with text to show (date format already provided (e.g. 'Dec 20, 2018 (6:00 UTC)')
          
IF dt_obj IS DATETIME OBJECT:
- date_format: str, format to display date (default: '%b %d, %Y (%H:%M UTC)')
    - example 1: '%b %d, %Y (%H:%M UTC)' could give 'Dec 20, 2018 (6:00 UTC)'
    - example 2: '%m-%d-%Y' could give '12-20-2018'
    
- method: method to place the date label (either 'anchor' for AnchoredText or 'manual' to place manually).
        (default: 'anchor')
- boxstyle: anchor box shape style (default: "round,pad=0.,rounding_size=0.2")
- facecolor: color of bounding box (default: 'black')
- edgecolor: color of bounding box edge (default: 'black')
- zorder: drawing order of date layer (default: 10)

IF METHOD = 'anchor':
- anchor_loc: anchor text location (default: 4)
- anchor_prop: anchor properties dictionary (default: {'size': 20, 'color':'white'})

IF METHOD = 'manual':
- x: x-location of figure extent to place date
- y: y-location of figure extent to place date
- textcolor: color oftext (default: 'white')
- fontsize: fontsize of text (defult: 15)

OUTPUT:
- input plot with added date label

DEPENDENCIES:
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
from datetime import datetime
from matplotlib.offsetbox import AnchoredText

Latest recorded update:
06-03-2022
    """

    
    assert method in ['anchor', 'manual'], f">>> method should be 'manual' or 'anchor', given: '{method}'"
    
    assert str(type(dt_obj)) in ["<class 'datetime.datetime'>", "<class 'str'>"], f">>> dt_obj should be datetime object or string, given: {str(type(dt_obj))}"
    
    
    # if given as datetime object, convert to specified date format
    if str(type(dt_obj)) == "<class 'datetime.datetime'>":
        date_text = dt_obj.strftime(date_format)
    
    # else, set date directly to given string object
    else:
        date_text = dt_obj
    
    
    # add text
    #---------
    if str(method) == 'anchor':
        at = AnchoredText(date_text, loc=anchor_loc, prop=anchor_prop)
        at.patch.set_boxstyle(boxstyle)
        at.patch.set_facecolor(facecolor)
        at.patch.set_edgecolor(edgecolor)
        at.zorder = zorder
        ax.add_artist(at)
    
    elif str(method) == 'manual':
        ax.text(x, y, date_text, 
                bbox=dict(boxstyle = boxstyle, facecolor=facecolor, edgecolor = edgecolor), 
                transform=ax.transAxes, fontsize=fontsize, 
                c=textcolor, verticalalignment='top', zorder = zorder);

        
        
#///////////////
#  add_mslp ///
#/////////////
#---------------------------------------------------------------------
# Add mslp contours to plot.
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
#---------------------------------------------------------------------
def add_mslp(fig, ax, Lons, Lats, msl, 
             
             levels = np.arange(980,1050,4), color_range = [980,1050], linewidths = 2, cmap = 'viridis', 
             
             fontsize = 15, cbar_label = 'Sea Level Pressure (hPa)',
             cbar_bottom = 0.15, cbar_top = 0.95,
             
             zorder = 5, ): 
    
    
    """Add mslp contours to plot.
    
INPUT:
- fig: cartopy figure
- ax: cartopy figure axis
- Lons: M x N lon grid for mslp data
- Lats: M x N lat grid for mslp data
- msl: M x N grid of mslp data
- levels: msl levels to plot
- color_range: Mx1 (M=2,3,4) lists of floats for color scale normalization, either as:
-              --> [vmin, vmax] (default: [980, 1050])
-              --> [vmin, midpoint, vmax] 
-              --> [vmin, midpoint1, midpoint2, vmax]
- linewidths: contour linewidths 
- cmap: colormap (default: 'viridis')
- fontsize: tick label and colorbar label fontsize (default: 15)
- cbar_label: colorbar label (default: 'Sea Level Pressure (hPa)')
- cbar_bottom: bottom height of cbar (as fraction of total height, default: 0.15)
- cbar_top: top height of cbar (as fraction of total height, default: 0.95)
- zorder: drawing order of layer (default: 5)

OUTPUT:
- input plot with added date label

DEPENDENCIES:
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors

Latest recorded update:
06-03-2022
    """
    
    
    
    # determine color normalization
    #------------------------------
    class TwopointNormalize(matplotlib.colors.Normalize):
        def __init__(self, vmin=None, vmax=None, vmid1=None, vmid2=None, clip=False):
            self.vmid1 = vmid1
            self.vmid2 = vmid2
            super().__init__(vmin, vmax, clip)

        def __call__(self, value, clip=None):
            # I'm ignoring masked values and all kinds of edge cases to make a
            # simple example...
            x, y = [self.vmin, self.vmid1, self.vmid2, self.vmax], [0, 0.33,0.66, 1]
            return np.ma.masked_array(np.interp(value, x, y))
    if len(color_range)==2:
        norm=matplotlib.colors.Normalize(vmin=color_range[0], vmax=color_range[1])
    elif len(color_range)==3:
        norm=matplotlib.colors.TwoSlopeNorm(vmin=color_range[0], vcenter=color_range[1],vmax=color_range[2])
    else:
        norm = TwopointNormalize(vmin  = color_range[0], vmid1 = color_range[1], 
                                 vmid2 = color_range[2], vmax  = color_range[3])
        
    
    # add MSLP contours
    #------------------
    Pcont = ax.contour(Lons, Lats, msl, levels=levels, 
                       cmap=cmap, norm = norm, linewidths = linewidths, 
                       transform=ccrs.PlateCarree(), zorder=zorder)
    
    # colorbar parameters
    #-------------------
    # normalize colorbar scale by provide min/max values    
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    # add colorbar
    #-------------
    # create ticks for 5 mslp levels
    y_bottom = ax.get_position().y0
    y_height = ax.get_position().height
    
    cbar = plt.colorbar(sm, 
                        ticks=np.around(np.linspace(color_range[0],color_range[-1],5),4),
                        cax=fig.add_axes([ax.get_position().x1+0.01, 
                                          y_bottom + cbar_bottom*y_height,
                                          0.02, 
                                          (cbar_top-cbar_bottom)*y_height]))
  
    # add colorbar with label
    #------------------------
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.ax.set_ylabel(cbar_label, fontsize = fontsize);

    
# NEW VERSION BELOW
# #//////////////////
# #  add_vectors ///
# #////////////////
# #---------------------------------------------------------------------
# # Add vectors to figure (automatically fix cartopy vector angle issue).
# #---------------------------------------------------------------------
# # DEPENDENCIES
# import numpy as np
# import matplotlib as mpl
# from matplotlib import pyplot as plt
# import matplotlib.colors
# import cartopy
# import cartopy.crs as ccrs
# #---------------------------------------------------------------------
# def add_vectors(fig, ax, Lons, Lats, u, v, regrid = [1,1], 
#               color = (0, 0, 0), scale=150, width = 0.005, 
#               headwidth = 3, headaxislength = 4, headlength = 4, minshaft = 1, minlength = 1, 
#               linewidth = 0, alpha=1, angles = 'uv', pivot = 'mid',
#               zorder_winds = 7, zorder_key = 6, quiv_key = (10, [1.1, 0.035],[1.02, 0.052], '            10 m/s', 15, 'k')): 
    
#     """Add vectors to figure (automatically fix cartopy vector angle issue).
    
# INPUT:
# - fig: cartopy figure
# - ax: cartopy figure axis
# - Lons: M x N lon grid for vectors
# - Lats: M x N lat grid for vectors

# THESE WILL BE RESCALED USING function: fix_cartopy_vectors
# - u: M x N grid for horizontal components of vectors
# - v: M x N grid for vertical components of vectors

# - regrid: regridding density, either specified as
#          --> list of length 1 specifying density to use cartopy regridding (e.g. [12])
#          --> list of length 2 specifying spacing along each grid direction (e.g. [5,10])
#          --> default: [1,1] (no regridding)
#                  if plotting a zonal transect of vectors (only one latitude, set wDlat = None) 
# - color: color of arrow vectors (default: (0.9, 0.4, 0))
# - alpha: opacity of vectors (default: 1)
# - scale: Number of data units per arrow length unit (default: 150) 
# - width: Shaft width in arrow units (default: 0.005)
# - headwidth: Head width as multiple of shaft width (default: 3)
# - headaxislength: Head length at shaft intersection (default: 4)
# - headlengt: Head length as multiple of shaft width (default:  4) 
# - minshaft: Length below which arrow scales, in units of head length (default: 1)
# - minlength: Minimum length as multiple of shaft width; if less than this, plot a dot of this diameter (default: 1)
# - linewidth: linewidth of arrow outline (default: 0)
# - angles: Method for determining the angle of the arrows (default: 'uv')
# - pivot: The part of the arrow that is anchored to the X, Y grid (default:  'mid')
# - zorder_winds: drawing order of winds layer (default: 7)
# - zorder_key: drawing order of wind key layer (default: 6)
# - quiv_key: key arrow size, text, fontsize, fontcolor, and position of key arrow and text box as 
#             (arrowsize, [Xarrow, Yarrow],[Xtext, Ytext], text, fontsize, fontcolor)
#             or set = None if no key is desired 
#             (default: (10, [1.1, 0.035],[1.02, 0.052], '            10 m/s', 15, 'k'))

# OUTPUT:
# - input plot with added date label

# DEPENDENCIES:
# import numpy as np
# import matplotlib as mpl
# from matplotlib import pyplot as plt
# import matplotlib.colors
# import cartopy
# import cartopy.crs as ccrs

# Latest recorded update:
# 07-14-2022
#     """
    
#     # plot vectors
#     #-------------
#     # set vector density spacing with regrid
#     # set arrow color as specified
    
#     # else plot with 2D grid spacing
#     #-------------------------------
#     if len(regrid) == 2: 
#         winds = ax.quiver(Lons[::regrid[0],::regrid[1]], Lats[::regrid[0],::regrid[1]], 
#                           *fix_cartopy_vectors(u[::regrid[0],::regrid[1]], v[::regrid[0],::regrid[1]], Lats[::regrid[0],::regrid[1]]), 
#                           transform=ccrs.PlateCarree(), 
#                           color = color, angles = angles, 
#                           scale=scale, width = width, headwidth = headwidth, 
#                           headaxislength = headaxislength, headlength = headlength,
#                           linewidth = linewidth, alpha=alpha,
#                           minshaft = minshaft, minlength = minlength, pivot = pivot, zorder = zorder_winds)

#     else:
#         winds = ax.quiver(Lons, Lats, *fix_cartopy_vectors(u, v, Lats), 
#                           transform=ccrs.PlateCarree(), regrid_shape=regrid[0],
#                           color = color, angles = angles, 
#                           scale=scale, width = width, headwidth = headwidth, 
#                           headaxislength = headaxislength, headlength = headlength, 
#                           linewidth = linewidth, alpha=alpha,
#                           minshaft = minshaft, minlength = minlength, pivot = pivot, zorder = zorder_winds)
    
    
#     # add vector key
#     #---------------
#     if quiv_key != None:
#         ax.quiverkey(winds, X=quiv_key[1][0], Y=quiv_key[1][1], U=quiv_key[0], label = '', labelpos='E')
#         ax.text(quiv_key[2][0], quiv_key[2][1], quiv_key[3], transform=ax.transAxes, fontsize=quiv_key[4], color=quiv_key[5],
#             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=1), zorder = zorder_key)
        
        
        
#//////////////////
#  add_vectors ///
#////////////////
#---------------------------------------------------------------------
# Add vectors to figure (automatically fix cartopy vector angle issue).
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
import cartopy
import cartopy.crs as ccrs
from LIB_geo_plot import fix_cartopy_vectors
#---------------------------------------------------------------------
def add_vectors(fig, ax, Lons, Lats, u, v, regrid = [1,1], 
              color = (0, 0, 0), scale=150, width = 0.005, 
              headwidth = 3, headaxislength = 4, headlength = 4, minshaft = 1, minlength = 1, 
              linewidth = 0, alpha=1, angles = 'uv', pivot = 'mid',
              zorder = 7, zorder_key = 6, quiv_key = (10, [1.1, 0.035],[1.02, 0.052], '            10 m/s', 15, 'None', 'None')): 
    
    """Add vectors to figure (automatically fix cartopy vector angle issue).
    
INPUT:
- fig: cartopy figure
- ax: cartopy figure axis
- Lons: M x N lon grid for vectors
- Lats: M x N lat grid for vectors

THESE WILL BE RESCALED USING function: fix_cartopy_vectors
- u: M x N grid for horizontal components of vectors
- v: M x N grid for vertical components of vectors

- regrid: regridding density, either specified as
         --> list of length 1 specifying density to use cartopy regridding (e.g. [12])
         --> list of length 2 specifying spacing along each grid direction (e.g. [5,10])
         --> default: [1,1] (no regridding)
                 if plotting a zonal transect of vectors (only one latitude, set wDlat = None) 
- color: color of arrow vectors (default: (0.9, 0.4, 0))
- alpha: opacity of vectors (default: 1)
- scale: Number of data units per arrow length unit (default: 150) 
- width: Shaft width in arrow units (default: 0.005)
- headwidth: Head width as multiple of shaft width (default: 3)
- headaxislength: Head length at shaft intersection (default: 4)
- headlengt: Head length as multiple of shaft width (default:  4) 
- minshaft: Length below which arrow scales, in units of head length (default: 1)
- minlength: Minimum length as multiple of shaft width; if less than this, plot a dot of this diameter (default: 1)
- linewidth: linewidth of arrow outline (default: 0)
- angles: Method for determining the angle of the arrows (default: 'uv')
- pivot: The part of the arrow that is anchored to the X, Y grid (default:  'mid')
- zorder: drawing order of vector layer (default: 7)
- zorder_key: drawing order of wind key layer (default: 6)
- quiv_key: key arrow size, position of key arrow and text box, text, fontsize, textbox edgecolor, textbox facecolor, as 
            (arrowsize, [Xarrow, Yarrow],[Xtext, Ytext], text, fontsize, box_edgecolor, box_facecolor)
            or set = None if no key is desired 
            (default: (10, [1.1, 0.035],[1.02, 0.052], '            10 m/s', 15, 'None', 'None'))

OUTPUT:
- input plot with added date label

DEPENDENCIES:
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
import cartopy
import cartopy.crs as ccrs

Latest recorded update:
09-21-2022
    """
    
    # plot vectors
    #-------------
    # set vector density spacing with regrid
    # set arrow color as specified
    
    # else plot with 2D grid spacing
    #-------------------------------
    if len(regrid) == 2: 
        winds = ax.quiver(Lons[::regrid[0],::regrid[1]], Lats[::regrid[0],::regrid[1]], 
                          *fix_cartopy_vectors(u[::regrid[0],::regrid[1]], v[::regrid[0],::regrid[1]], Lats[::regrid[0],::regrid[1]]), 
                          transform=ccrs.PlateCarree(), 
                          color = color, angles = angles, 
                          scale=scale, width = width, headwidth = headwidth, 
                          headaxislength = headaxislength, headlength = headlength,
                          linewidth = linewidth, alpha=alpha,
                          minshaft = minshaft, minlength = minlength, pivot = pivot, zorder = zorder)

    else:
        winds = ax.quiver(Lons, Lats, *fix_cartopy_vectors(u, v, Lats), 
                          transform=ccrs.PlateCarree(), regrid_shape=regrid[0],
                          color = color, angles = angles, 
                          scale=scale, width = width, headwidth = headwidth, 
                          headaxislength = headaxislength, headlength = headlength, 
                          linewidth = linewidth, alpha=alpha,
                          minshaft = minshaft, minlength = minlength, pivot = pivot, zorder = zorder)
    
    
    # add vector key
    #---------------
    if quiv_key != None:
        key = ax.quiverkey(winds, X=quiv_key[1][0], Y=quiv_key[1][1], U=quiv_key[0], label = '', labelpos='E')
        key.set(zorder = zorder_key)
        ax.text(quiv_key[2][0], quiv_key[2][1], quiv_key[3], transform=ax.transAxes, fontsize=quiv_key[4],
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor=quiv_key[6], edgecolor=quiv_key[5], alpha=1), zorder = zorder_key+1)


#////////////////////////
#  plot_scalar_mesh  ///
#////////////////////// 
#---------------------------------------------------------------------
# function for plotting regular or geographic data with pcolormesh
#---------------------------------------------------------------------
# DEPENDENCIES:
import matplotlib.colors
import numpy as np, numpy.ma as ma
import cartopy, cartopy.crs as ccrs
import matplotlib.cm as cm
#---------------------------------------------------------------------

def plot_scalar_mesh(fig, ax, x_grid, y_grid, scalar_grid, geo_grid = True,
                     cmap = 'RdBu', cmap_norm = 'auto', cmap_type = 'continuous',
                     shading='nearest', zorder = 0):
    """Function for plotting regular or geographic data with pcolormesh. Many commands are automated in this function to save space in routine where it is used.

INPUT: 
- fig: figure to which scalar mesh will be added
- ax: figure axis to which scalar mesh will be added
- geo_grid: bool, whether or not grid is showing geographic data
    True: then x_grid is longitudes and y_grid is latitudes,
    data will be transformed using ccrs.PlateCarree() (default)
    False: data will not be transformed 
- x_grid: MxN grid of longitudes
- y_grid: MxN grid of latitudes
- scalar_grid: MxN grid of scalar data values
- cmap: colormap to use (default: 'RdBu')
- cmap_norm: normalization for data in colormap (default: 'auto')
    'auto' where normalization generated from scalar_grid data
    Mx1 lists of floats (M=2,3,4 when cmap_type == 'continuous'), 
    either as:[vmin, vmax], [vmin, midpoint, vmax], or [vmin, midpoint1, midpoint2, vmax]
    (M=any length when cmap_type == 'discrete')
- cmap_type: either 'discrete' or 'continuous' (default: 'continuous')
    'discrete' plots mesh as discretized data along boundaries given by cmap_norm
    'continuous' plots mesh as continuous colormap 
- shading: pcolormesh shading method (default: 'nearest')
- zorder: zorder of mesh layer


OUTPUT:
- fig, ax: figure and figure axis to which scalar mesh was added
- scalar_mesh: pcolormesh plot output

DEPENDENCIES:
import matplotlib.colors
import numpy as np, numpy.ma as ma
import cartopy, cartopy.crs as ccrs
import matplotlib.cm as cm
* also uses homemade TwopointNormalize class

Latest recorded update:
04-20-2022
    """

    # COLORMAP NORMALIZATION
    #=======================
    if str(cmap_type) == 'discrete':
        divnorm = matplotlib.colors.BoundaryNorm(cmap_norm, cmap.N)
    else:
        # automatically scale colormap 
        # based off of scalar_grid data
        if str(cmap_norm) == 'auto':
            divnorm = matplotlib.colors.TwoSlopeNorm(vmin=np.nanmin(scalar_grid), 
                                                     vcenter=(np.nanmax(scalar_grid)-np.nanmin(scalar_grid))/2+np.nanmin(scalar_grid), 
                                                     vmax=np.nanmax(scalar_grid))
        # create colormap normalization 
        # based on given input camp_norm values
        elif len(cmap_norm) == 2:
            divnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmap_norm[0], 
                                                     vcenter=(cmap_norm[1]+cmap_norm[0])/2, 
                                                     vmax=cmap_norm[1])
        elif len(cmap_norm) == 3:
            divnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmap_norm[0], 
                                                     vcenter=cmap_norm[1], 
                                                     vmax=cmap_norm[2])
        elif len(cmap_norm) == 4:
            divnorm = TwopointNormalize(vmin=cmap_norm[0], 
                                        vmid1=cmap_norm[1], 
                                        vmid2=cmap_norm[2], 
                                        vmax=cmap_norm[3])
        else:
            divnorm = cmap_norm

    # PLOT SCALAR MESH
    #=================
    # transform coordinates if geographic data
    if geo_grid == True:
        scalar_mesh=ax.pcolormesh(x_grid, y_grid, scalar_grid, 
                   cmap=cmap, norm=divnorm, transform=ccrs.PlateCarree(), shading=shading,zorder=zorder)
    else:
        scalar_mesh=ax.pcolormesh(x_grid, y_grid, scalar_grid, 
                   cmap=cmap, norm=divnorm, shading=shading, zorder=zorder)
    
    return fig, ax, scalar_mesh


#/////////////////////
#  plot_contours  ///
#///////////////////
#---------------------------------------------------------------------
# Function for plotting regular or geographic data contours
#---------------------------------------------------------------------
# DEPENDENCIES:
import matplotlib.colors
import numpy as np, numpy.ma as ma
import cartopy, cartopy.crs as ccrs
import matplotlib.cm as cm
#---------------------------------------------------------------------

def plot_contours(fig, ax, x_grid, y_grid, scalar_grid, geo_proj = 'None',
                  levels = 'auto', lw = 1, color = 'black',  
                  label_contours = False, manual_labels = False,
                  x_labels = np.array([180, 180, 180, 180, 180]), 
                  y_labels = np.array([65, 70, 75, 80, 85]), 
                  hidecont_belowlabel= True,
                  pixelspace_aroundlabel=30, labelsize=10, zorder = 1):
                      
    """Function for plotting regular or geographic data contours. Many commands are automated in this function to save space in routine where it is used.

INPUT: 
- fig: figure to which scalar mesh will be added
- ax: figure axis to which scalar mesh will be added
- x_grid: MxN grid of longitudes
- y_grid: MxN grid of latitudes
- scalar_grid: MxN grid of scalar data values
- geo_proj: either 'None' (default) or set to cartopy map projection (e.g. ccrs.NorthPolarStereo(central_longitude=203.5)) in which case x_grid is longitudes and y_grid is latitudes, and data will be transformed using ccrs.PlateCarree() 
- levels: contour levels to plot or 'auto' (default: 'auto', automatically generate from grid)
- lw: contour linewidth (default: 1)
- color: contour color (default: 'black')
- label_contours: bool, whether or not to label contours (default: False)
- manual_labels: bool, whether or not to add manual labels to contours (default: False)
- hidecont_belowlabel: bool, whether or not to remove contour below labels (default: True)
- pixelspace_aroundlabel: pixel spacing around label (default: 30) 
- labelsize: label size (default: 10)
- zorder: zorder of mesh layer
if label_contours == True and manual_labels == True
- x_labels: array of x or longitude positions to add to contour labels (default along dateline)
- y_labels: array of y or latitude positions to add to contour labels (default values 65-80)

OUTPUT:
- fig, ax: figure and figure axis to which contour plot was added

DEPENDENCIES:
import matplotlib.colors
import numpy as np, numpy.ma as ma
import cartopy, cartopy.crs as ccrs
import matplotlib.cm as cm

Latest recorded update:
04-20-2022
    """
    
    # make plot
    #=======================    
    # if not specified, plot data without projection
    if str(geo_proj) == 'None':
        if str(levels) != 'auto':
            CS =  ax.contour(x_grid, y_grid, scalar_grid, levels = levels, 
                         colors=color,linewidths=cont_lw, zorder=zorder)
        else:
            CS =  ax.contour(x_grid, y_grid, scalar_grid,
                         colors=color,linewidths=cont_lw, zorder=zorder)
    # if specified, project data onto map    
    else:
        if str(levels) != 'auto':
            CS =  ax.contour(x_grid, y_grid, scalar_grid, levels = levels, 
                             colors=color,linewidths=lw, transform = ccrs.PlateCarree(), zorder=zorder)
        else:
            CS =  ax.contour(x_grid, y_grid, scalar_grid,
                             colors=color,linewidths=lw, transform = ccrs.PlateCarree(), zorder=zorder)
            
    # add labels to contours
    #=======================
    if label_contours == True:
        # check whether to add manual labels
        if manual_labels == True:
            # if not specified, create manual points without transforming
            if geo_proj == False:
                manual_points = []
                for ii in range(len(x_labels)):
                    manual_points.append((x_labels[ii],y_labels[ii]))
            # if specified, transform provided manual points
            else:
                Coords = geo_proj.transform_points(ccrs.PlateCarree(), x_labels, y_labels)
                manual_points=[]
                for spot in Coords[:,0:2]:
                    manual_points.append((spot[0],spot[1]))

            # label contours with manual labels 
            ax.clabel(CS, CS.levels, inline=hidecont_belowlabel, 
                      inline_spacing=pixelspace_aroundlabel, 
                      fontsize=labelsize, manual=manual_points)
        else:
            # label contours without manual labels 
            ax.clabel(CS, CS.levels, inline=hidecont_belowlabel, 
                      inline_spacing=pixelspace_aroundlabel, 
                      fontsize=labelsize)
                  
    return fig, ax




#//////////////////
#  add_LFice   ///
#////////////////
#---------------------------------------------------------------------
# Plot 'landfast ice' perimeter around Alaskan coast from bathymetry data
#---------------------------------------------------------------------
# DEPENDENCIES:
import xarray as xr
import numpy as np, numpy.ma as ma
from matplotlib import pyplot as plt
import cartopy, cartopy.crs as ccrs
#---------------------------------------------------------------------


def add_LFice(ax, lat_range = [67, 71.5], lon_range = [190, 233], spacing = 5, depth_cutoff = -20,
              zorder = 8, color = [0.4,0.4,0.4], size = 0.01, marker = '.', bath_data = '/Volumes/Jewell_EasyStore/AcrticBathymetry_IBCAO_2022/gebco_2022_n90.0_s60.0_w0.0_e360.0.nc'):

    """
    Add 'landfast ice' cover perimeter around Alaskan coast to plot, using bathymetry data \n
    20 m depth and shallower based on Mahoney et a. (2007), doi: 10.1029/2006JC003559
    
    
    INPUT:
    - ax: figure axis
    - bath_data: Beaufort Sea bathymetry data (default: '/Volumes/Jewell_EasyStore/AcrticBathymetry_IBCAO_2022/gebco_2022_n90.0_s60.0_w0.0_e360.0.nc')
    - lat_range: range of lats to crop as [min, max] (default:[67, 71.5])
    - lon_range: range of lons to crop as [min, max] (default:[190, 233])
    - spacing: step size in either direction used to grab data, value of 1 selects all data points (default: 5)
    - depth_cutoff: deepest bathymetry (m) to include as part of LF ice (default: -20)
    - color: color of LF ice region (default: [0.4,0.4,0.4])
    - point_size: size of scatter points to plot (default: 0.01)
    - zorder: plot layer (default: 8)
    - marker: marker type of scatter plot (default: '.')
    
    DEPENDENCIES:
    import xarray as xr
    import numpy as np, numpy.ma as ma
    from matplotlib import pyplot as plt
    import cartopy, cartopy.crs as ccrss

    Latest recorded update:
    12-14-2022
    """
    
    # import elevation data
    ds = xr.open_dataset(bath_data)
    ds.close()
    ds_crop = ds.sel(lat=slice(lat_range[0],lat_range[1]),lon=slice(lon_range[0],lon_range[1]))
    
    # grab geo coordinates and elevation data
    lat_elev = ds_crop.lat
    lon_elev = ds_crop.lon
    elev = ds_crop.elevation
    
    # apply reduction in point density if specified
    lon, lat = np.meshgrid(lon_elev[::spacing], lat_elev[::spacing])
    elev = elev[::spacing, ::spacing]
    lon_shallower = lon[elev>=depth_cutoff]
    lat_shallower = lat[elev>=depth_cutoff]
    
    ax.scatter(lon_shallower,lat_shallower, s=size, color=color, marker=marker, transform = ccrs.PlateCarree(), zorder=zorder)
   
