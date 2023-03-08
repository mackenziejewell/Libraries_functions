#/////////////////////////
#  TwopointNormalize  ///
#///////////////////////
#---------------------------------------------------------------------
# class for normalizing colormap based off two midpoints
#---------------------------------------------------------------------
# DEPENDENCIES:
import matplotlib.colors
import numpy.ma as ma
#---------------------------------------------------------------------

class TwopointNormalize(matplotlib.colors.Normalize):
    
    """Class for normalizing colormap based off two midpoints.

INPUT: 
- vmin: min value
- vmid1: lower midpoint value
- vmid2: higher midpoint value
- vmax: max value

OUTPUT:
- normalization scaling [vmin, vmid1, vmid2, vmax] to [0, 1/3, 2/3, 1] of colormap

DEPENDENCIES:
import matplotlib.colors
import numpy.ma as ma

Latest recorded update:
04-20-2022
    """
        
    def __init__(self, vmin=None, vmax=None, vmid1=None, vmid2=None, clip=False):
        self.vmid1 = vmid1
        self.vmid2 = vmid2
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.vmid1, self.vmid2, self.vmax], [0, 0.33,0.66, 1]
        return np.ma.masked_array(np.interp(value, x, y))
    
       
    
#////////////////////
#  add_colorbar  ///
#//////////////////
#---------------------------------------------------------------------
# add colorbar to plot, with functions to reduce clutter
#---------------------------------------------------------------------
# DEPENDENCIES:
import matplotlib.colors
import numpy as np, numpy.ma as ma
import cartopy, cartopy.crs as ccrs
import matplotlib.cm as cm
#---------------------------------------------------------------------

def add_colorbar(fig, ax, colorbar_input, cb_placement = 'left', cb_orientation = 'auto', 
                 cb_width = 'auto',  cb_length_fraction = [0,1], cb_pad = 0, 
                 cb_ticks = 'auto', cb_ticklabels = 'auto', cb_extend='neither', cb_label=' ', 
                 cb_labelsize = 12, draw_edges=False, edge_params=['k',2], suppress_prints=True):
    
    """Function for plotting colorbar along edge of figure axis.

INPUT: 
- fig: figure to which colorbar will be added
- ax: figure axis to which colorbar will be added
- colorbar_input: either specify [matplotlib.collections.QuadMesh], pmatplotlib.cm.ScalarMappable] (from pcolormesh plot output),
                    [cartopy.mpl.contour.GeoContourSet] (from countourf output),
                  or specify [cmap, norm] 
                   where cmap is matplotlib cmap (e.g. 'RdBu')
                   where norm is matplotlib.colors normlalization instance (e.g. made from TwoSlopeNorm)
- cb_placement: location/orientation of colorbar, as 'left' (default),'right','top','bottom'
- cb_orientation: orientation ('horizontal' or 'vertical') of colorbar. Set to 'auto' (default) to 
                  pick automatically based off its cb_placement
- cb_width: colorbar width (default: 'auto', which makes it 1/20 figure width)
- cb_length_fraction: beginning and end position of colorbar along axis as [begin, end], as fraction of axis length 
                      (default: [0,1] for cbar to stretch along full axis)
- cb_pad: pad between plot and colorbar (default: 0)
- cb_ticks: colorbar ticks. 'auto' (default) selects automatically from data, or provide ticks as list (e.g. [1,2,3])
- cb_ticklabels:  colorbar tick labels
             'auto' (default) selects automatically from data, or provide ticks as list (e.g. ['<1','2','>3'])
              if providing list, must match number of provided cb_ticks
- cb_extend: end cap style for colorbar (to address out-of-range values), either:
           --> 'neither': (default) flat ends at either end
           --> 'min': arrow at min end of colorbar
           --> 'max': arrow at max end of colorbar
           --> 'both': arrow at both ends of colorbar
- cb_label: colorbar label (string), default is empty string
- cb_labelsize: colorbar label and tick fontsize
- draw_edges: bool, whether or not to draw outline around colorbar (default:False)
- edge_params: color and linewidth for cbar edges if drawn, as [edgecolor, edgelinewidth] (default: ['k',2])


OUTPUT:
- fig: figure to which colorbar was added
- ax: figure axis to which colorbar was added

DEPENDENCIES:
import matplotlib.colors
import numpy as np, numpy.ma as ma
import cartopy, cartopy.crs as ccrs
import matplotlib.cm as cm

Latest recorded update:
10-10-2022
    """
    
    
    # determine type of colorbar input 
    #=================================
    # if colorbar_input is [QuadMesh] from plot output
    if len(colorbar_input) == 1:
        if isinstance(colorbar_input[0], matplotlib.collections.QuadMesh) or isinstance(colorbar_input[0], cartopy.mpl.contour.GeoContourSet) or isinstance(colorbar_input[0], matplotlib.cm.ScalarMappable):
            
            CB_INPUT = colorbar_input[0]
        else:
            if suppress_prints==False:
                print('colorbar_input is not type matplotlib.collections.QuadMesh nor cartopy.mpl.contour.GeoContourSet')
    # if colorbar_input is [cmap, norm]
    elif len(colorbar_input) == 2:
        CB_INPUT = cm.ScalarMappable(norm=colorbar_input[1], cmap=colorbar_input[0])
    else:
        if suppress_prints==False:
            print('unrecognized colorbar_input, should be of length 1 or 2')
    
    # generate plot axes
    #=================================
    # get plot axes corner coordinates
    plot_axis_coords = ax.get_position().get_points()
    ax_x0 = plot_axis_coords[0][0]
    ax_x1 = plot_axis_coords[1][0]
    ax_y0 = plot_axis_coords[0][1]
    ax_y1 = plot_axis_coords[1][1]
    
    # grab desored fractional lengths of colorbar
    #============================================
    cb_L_i = cb_length_fraction[0]
    cb_L_f = cb_length_fraction[1] 

    # set widths of colorbar based of specification or 1/10 figure width
    if str(cb_width) == 'auto':
        if str(cb_placement) == 'top' or str(cb_placement) == 'bottom':
            WIDTH = 0.05*(ax_y1-ax_y0)
        else:
            WIDTH = 0.05*(ax_x1-ax_x0)
    else:
        WIDTH = cb_width

    # generate colorbar axis based off desired edge of placement
    if str(cb_placement) == 'left':  
        cbar_ax = fig.add_axes([ax_x0-(WIDTH+cb_pad), ax_y0+(cb_L_i*(ax_y1-ax_y0)), WIDTH, (ax_y1-ax_y0)*(cb_L_f-cb_L_i)])
    elif str(cb_placement) == 'right':
        cbar_ax = fig.add_axes([ax_x1+cb_pad, ax_y0+(cb_L_i*(ax_y1-ax_y0)), WIDTH, (ax_y1-ax_y0)*(cb_L_f-cb_L_i)])
    elif str(cb_placement) == 'top':
        cbar_ax = fig.add_axes([ax_x0+(cb_L_i*(ax_x1-ax_x0)), ax_y1+cb_pad, (ax_x1-ax_x0)*(cb_L_f-cb_L_i), WIDTH])
    else:
        cbar_ax = fig.add_axes([ax_x0+(cb_L_i*(ax_x1-ax_x0)), ax_y0-(WIDTH+cb_pad), (ax_x1-ax_x0)*(cb_L_f-cb_L_i), WIDTH])
        
    # set colorbar orientation from its placement
    if str(cb_orientation) == 'auto':
        if str(cb_placement) == 'top' or str(cb_placement) == 'bottom':
            cb_orientation = 'horizontal'
        else:
            cb_orientation = 'vertical'

    # make colorbar and place labels
    #====================================
    # if colorbar ticks not provided, automatically place ticks
    if str(cb_ticks) == 'auto':
        cbar = fig.colorbar(CB_INPUT,cax=cbar_ax, 
                            orientation=cb_orientation, extend=cb_extend, drawedges=draw_edges)

    # if colorbar ticks not provided, place as specified
    else:
        cbar = fig.colorbar(CB_INPUT,cax=cbar_ax, 
                            orientation=cb_orientation, extend=cb_extend, ticks=cb_ticks, drawedges=draw_edges)
        # place tick labels if specified
        if str(cb_ticklabels)!='auto':
            if str(cb_placement) == 'top' or str(cb_placement) == 'bottom':
                cbar.ax.set_xticklabels(cb_ticklabels) 
            else:
                cbar.ax.set_yticklabels(cb_ticklabels) 

    # if including edge border around cbar, specify its linewidth and color            
    if draw_edges==True:
        cbar.outline.set_color(edge_params[0])
        cbar.outline.set_linewidth(edge_params[1])
    
    # remove gray facecolor behind colorbar arrow cap if extending at either end
    if str(cb_extend) != 'neither':
        cbar.ax.set_facecolor('none')
        
    # set label and colorbar fontsize
    cbar.ax.tick_params(labelsize=cb_labelsize)
    cbar.set_label(cb_label, fontsize=cb_labelsize)

    # place ticks and label on correct side of colorbar for various colorbar orientations
    if str(cb_placement) == 'top' or str(cb_placement) == 'bottom':
        cbar_ax.xaxis.set_ticks_position(cb_placement)
        cbar_ax.xaxis.set_label_position(cb_placement)
    else:
        cbar_ax.yaxis.set_ticks_position(cb_placement)
        cbar_ax.yaxis.set_label_position(cb_placement)


    return fig, ax