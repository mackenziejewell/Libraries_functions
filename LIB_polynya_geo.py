# DEPENDENCIES:
from skimage.segmentation import find_boundaries
from scipy.ndimage import binary_fill_holes
import numpy as np, numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
from skimage.morphology import flood
from skimage.morphology import disk
import math

def save_to_grid(final_groups, polynyas,  keys, cellarea, xx, yy, grid, ice, 
                 land_no_holes, minimum_area = 0, show_plots = True):
    
    """Save final grouped polynyas to grid.

INPUT: 
- final_groups: updated dictionary of polynya groups.
- polynyas: dictinary of grouped polynya keys, groups increase from 1.
- keys: (M x N) grid of distinct keys at each coordinate
- cellarea: (M x N) grid of pixel areas
- xx: (N x 1) array of x coordinates
- yy: (M x 1) array of y coordinates
- grid: (M x N) array of categories
- ice: (M x N) boolean grid of ice pixels
- land_no_holes: (M x N) boolean grid, true over land areas.
- minimum_area: minimum_area (sq. km) requirement for saving final polynyas (default: 0)
- show_plots

OUTPUT:
- final_grid: final (M x N) grid, -2 = land, -1 = ice, 0 = open ocean, 1+ = polynyas
- large_polynyas: dictionary of grouped polynyas

DEPENDENCIES:
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from skimage.morphology import disk
from skimage.segmentation import find_boundaries
import math

homemade function: grab_footprint

Latest recorded update:
05-22-2024
    """
    

    final_grid = np.zeros_like(grid, dtype=int)
    final_grid[ice] = -1
    final_grid[land_no_holes.astype(int)==1] = -2

    large_polynyas = {}

    ppp = 1
    for group in final_groups.keys():

        # grab numbers of all polynyas in group
        all_p_keys = final_groups[group]

        # find all keys for all polynyas in group
        all_keys = np.array([], dtype=int)
        for p_key in all_p_keys:
            all_keys = np.append(all_keys, polynyas[p_key]['edge'])
            all_keys = np.append(all_keys, polynyas[p_key]['inner'])

        # calculate total polynya area
        total_polynya_area = 0
        for key in all_keys:
            total_polynya_area += cellarea[np.where(keys == key)][0]/(1000*1000)

        if total_polynya_area > minimum_area:

            large_polynyas[ppp] = {}
            large_polynyas[ppp]['keys'] = all_keys

            all_x = np.array([])
            all_y = np.array([])

            for key in all_keys:
                (ii, jj) = np.where(keys == key)
                ii, jj = ii[0], jj[0]

                all_y = np.append(all_y, yy[ii])
                all_x = np.append(all_x, xx[jj])

                final_grid[ii, jj] = ppp

            large_polynyas[ppp]['X'] = np.mean(all_x)
            large_polynyas[ppp]['Y'] = np.mean(all_y)

            ppp+=1

    if show_plots:
        # plotting
        #---------

        # custom colormap
        cmap = plt.cm.Paired              # get a specific colormap
        cmaplist = cmap.colors[:-2]                      # extract all colors
        X = math.ceil(len(large_polynyas.keys()) / len(cmaplist))
        cmaplistext = cmaplist
        for n in range(X-1):
            cmaplistext = cmaplistext + cmaplist
        customMap = matplotlib.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplistext, len(large_polynyas.keys()))

        fig, ax = plt.subplots(figsize=(4,6))
        ax.pcolormesh(xx, yy, final_grid, cmap=customMap)
        ax.pcolormesh(xx, yy, final_grid==0, cmap=matplotlib.colors.ListedColormap(['None', 'lightgray']))
        ax.pcolormesh(xx, yy, final_grid==-1, cmap=matplotlib.colors.ListedColormap(['None', 'white']))
        ax.pcolormesh(xx, yy, final_grid==-2, cmap=matplotlib.colors.ListedColormap(['None', 'gray']))

        for key in large_polynyas.keys():
            ax.text(large_polynyas[key]['X'], large_polynyas[key]['Y'], key, ha='center', va = 'center', weight='normal')
    
    
    return final_grid, large_polynyas



def merge_polynyas(keys, polynyas, polynya_no_marginal, radius = 4, quiet = False):

    """Merge polynyas within radius of one another.

INPUT: 
- keys: (M x N) grid of distinct keys at each coordinate
- polynyas: dictinary of grouped polynya keys, groups increase from 1.
- polynya_no_marginal: (M x N) boolean grid of polynyas not within marginal ice zone
- radius: radius to use in footprint (radius 1 has length of xx, yy, coordinate spacing).
- quiet: bool, whether or not to print statements.

OUTPUT:
- final_groups: updated dictionary of polynya groups.

DEPENDENCIES:
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from skimage.morphology import disk
from skimage.segmentation import find_boundaries

homemade function: grab_footprint

Latest recorded update:
05-22-2024
    """
    
    # 1: Loop through all polynyas and 
    # find any within radius of one another
    #---------------------------------------
    # create 50 km radius footprint to check
    footprint = disk(radius).astype(bool)
    ii_size = int((footprint.shape[0]-1)/2)
    jj_size = int((footprint.shape[1]-1)/2)

    polynyas_to_check = np.unique(polynya_no_marginal[polynya_no_marginal>0])
    grouped_polynyas = {}
    fully_checked_poly = np.array([], dtype=int)

    if not quiet:
        print(f'Of {len(polynyas.keys())} original polynyas, proceed with {len(polynyas_to_check)} non-marginal polynyas')

    for p_num in polynyas_to_check:

        # if it hasn't already been fully checked
        if p_num not in fully_checked_poly:

            # start new group of polynyas to save
            grouped_polynyas[p_num] = np.array([p_num])

            # find all edges of current polynya
            p_edges = polynyas[p_num]['edge']

            # record that polynya will be fully checked
            fully_checked_poly = np.append(fully_checked_poly, p_num)

            # loop through polynya edges and find all others within 50 km radius
            for pol_key in p_edges:

                (ii, jj) = np.where(keys == pol_key)
                ii, jj = ii[0], jj[0]
                # find values in radius 4 around ii, jj
                foot = footprint
                nearby_vals = grab_footprint(foot, polynya_no_marginal, ii, jj)

                # find all numbers of other polynyas within radius and add to group
                nonzero_diff_vals = np.unique(nearby_vals[(nearby_vals!=0)&(nearby_vals!=p_num)])
                for val in nonzero_diff_vals:
                    if val not in grouped_polynyas[p_num]:
                        grouped_polynyas[p_num] = np.append(grouped_polynyas[p_num], val)


    # 2: Consolidate into groups
    #---------------------------------------
    final_groups = {}
    fully_checked_poly = np.array([], dtype=int)

    for polynya in list(grouped_polynyas.keys()):

        if polynya not in fully_checked_poly:

            fully_checked_poly = np.append(fully_checked_poly, polynya)

            # if just individual polynya, move over to next dict
            if len(grouped_polynyas[polynya]) == 1:
                final_groups[polynya] = grouped_polynyas[polynya]

            # if grouped polynyas, merge groups
            else:

                # start new polynya group if needed
                if polynya not in list(final_groups.keys()):
                    final_groups[polynya] = grouped_polynyas[polynya]

                # iteratively add all attached polynyas on until number
                # of polynyas in group stops growing
                size_change = 100
                ii = 0
                while size_change > 0:

                    polynya_size = len(final_groups[polynya])

                    # loop through all polynyas in group to check which others its linked to
                    for other_poly in final_groups[polynya]:
                        if other_poly not in fully_checked_poly:
                            # add on all polynyas linked to other_poly
                            final_groups[polynya] = np.append(final_groups[polynya], grouped_polynyas[other_poly])
                            # remove repeats
                            final_groups[polynya] = np.unique(final_groups[polynya])     
                            # record that other_poly has been addresses
                            fully_checked_poly = np.append(fully_checked_poly, other_poly)

                    size_change = len(final_groups[polynya]) - polynya_size
                    polynya_size = len(final_groups[polynya])

                    # stop if somehow loops more than 100 times
                    ii+=1
                    if ii > 100:
                        size_change = 0

    if not quiet:
        print(f'>>>> consolidate polynyas into {len(final_groups.keys())} groups')
    
    return final_groups



def remove_marginal_polynyas(grid, polynya_flooded, ocean_edge, ice, land_no_holes, keys, radius = 4, 
                             show_plots = {'xx': [], 'yy': []}):

    """Remove polynyas within radius of ice edge.

INPUT: 
- grid: (M x N) category grid: 0 = land, 1 = open ocean, 2 = polynyas, 3 = ice
- polynya_flooded: (M x N) grid of grouped polynyas, #s indicate polynya group. 0 = no polynya.
- ocean_edge: (M x N) boolean grid of ocean edge pixels
- ice: (M x N) boolean grid of ice pixels
- land_no_holes: (M x N) boolean grid, true over land areas.
- keys: (M x N) grid of distinct keys at each coordinate
- radius: radius to use in footprint (radius 1 has length of xx, yy, coordinate spacing).
- show_plots: None if you don't want to show plots, or dict of plotting parameters 
    {'xx': (MxN) x-coordinates, 'yy': (MxN) xy-coordinates, 'date': date of data}

OUTPUT:
- polynya_no_marginal: updated (M x N) polynya_flooded grid, polynyas in MIZ removed
- updated_grid: (M x N) updated version of grid with marginal polynyas replaced as open ocean
- new_just_polynyas: (M x N) boolean grid of polynyas not within MIZ
- final_polynya_edges: (M x N) updated boolean grid of polynya edges
- _ocean_edge_: (M x N) updated boolean grid of ocean edges, including those in MIZ

DEPENDENCIES:
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from skimage.morphology import disk
from skimage.segmentation import find_boundaries

homemade function: grab_footprint

Latest recorded update:
05-22-2024
    """
    
    # create 50 km radius footprint to check
    footprint = disk(radius).astype(bool)
#     ii_size = int((footprint.shape[0]-1)/2)
#     jj_size = int((footprint.shape[1]-1)/2)

    # list of polynyas to remove
    marginal_polynyas = np.array([], dtype=int)

    marginal_grid = np.zeros_like(polynya_flooded)
    polynya_no_marginal = np.copy(polynya_flooded)

    _ocean_edge_ = ocean_edge

    # diff will record # of polynyas recorded as marginal on each iteration 
    diff = 100
    num_remove = 0
    
    # stop running after 100 iterations if for some reason does not stop
    iteration = 0

    # continue looping while new marginal polynyas being found
    while diff > 0:
    
        if iteration == 100:
            break

        elif iteration == 0:
            # run footprint around each pixel in ocean edge
            ocean_edge_keys = keys[_ocean_edge_]

        else:
            # run footprint around each pixel in ocean edge
            ocean_edge_keys = keys[new_ocean_edge_==1]

        # loop through current ocean edge keys and check whether any
        # polynyas are within radius of pixel
        for ocean_key in ocean_edge_keys:

            # find values in radius around ii,jj 
            (ii, jj) = np.where(keys == ocean_key)
            ii, jj = ii[0], jj[0]
            foot = footprint
            nearby_vals = grab_footprint(foot, polynya_flooded, ii, jj)
            marginal_polynyas = np.append(marginal_polynyas, nearby_vals)

        # list every-non ocean polynya key only once
        marginal_polynyas = np.unique(marginal_polynyas[marginal_polynyas>0])

        # update marginal, non-marginal grids
        for p_key in marginal_polynyas:
            polynya_no_marginal[(polynya_flooded == p_key)] = 0
            marginal_grid[(polynya_flooded == p_key)] = 1

        # update grid and new ice edge
        updated_grid = grid - marginal_grid

        # find ocean edge along ice pack
        new_just_polynyas = (updated_grid==2)
        # polynya edges don't include boundaries with land this way, to be used for next_ocean_edge_ calc
        polynya_edges = (new_just_polynyas == 1) & ((find_boundaries(ice, mode='outer') - land_no_holes.astype(int))==1)
        next_ocean_edge_ = ((find_boundaries(ice, mode='outer') - land_no_holes.astype(int))==1 - polynya_edges)==1
        new_ocean_edge_ = next_ocean_edge_.astype(int) - _ocean_edge_.astype(int)
        _ocean_edge_ = next_ocean_edge_

        diff = len(marginal_polynyas) - num_remove
        num_remove = len(marginal_polynyas)

        iteration+=1
        
        
    final_polynya_edges = find_boundaries(new_just_polynyas, mode='inner')
    
    if show_plots != None:
        
        xx = show_plots['xx']
        yy = show_plots['yy']
        
        # plotting
        #---------
        fig, ax = plt.subplots(figsize=(4.5,6))
        bounds = np.array([-0.5, 0.5, 1.5, 2.5,3.5])
        grid_norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=len(bounds)-1)
        grid_camp = matplotlib.colors.ListedColormap(['darkgray','lightsteelblue','lightgray', 'white',])
        mesh = plt.pcolormesh(xx, yy, updated_grid, norm = grid_norm, cmap=grid_camp)
        ax.pcolormesh(xx, yy, _ocean_edge_, 
                       cmap=matplotlib.colors.ListedColormap(['None', 'skyblue']))
        ax.pcolormesh(xx, yy, ocean_edge, 
                       cmap=matplotlib.colors.ListedColormap(['None', 'blue']))
        ax.pcolormesh(xx, yy, ma.masked_where(polynya_no_marginal==0, polynya_no_marginal), cmap='cmo.amp')
        ax.pcolormesh(xx, yy, final_polynya_edges, 
                       cmap=matplotlib.colors.ListedColormap(['None', 'lightcoral']))


    return polynya_no_marginal, updated_grid, new_just_polynyas, final_polynya_edges, _ocean_edge_



def group_continuous_polynyas(polynya_edge, just_polynyas, keys, 
                              show_plots = {'xx': [], 'yy': [], 'grid': []}):
    
    """Group polynyas into connected features.

INPUT: 
- polynya_edge: (M x N) boolean grid of polynya edge pixels
- just_polynyas: (M x N) boolean grid of polynya pixels
- keys: (M x N) grid of distinct keys at each coordinate
- show_plots: None if you don't want to show plots, or dict of plotting parameters 
    {'xx': (MxN) x-coordinates, 'yy': (MxN) xy-coordinates, 'date': date of data}

OUTPUT:
- polynyas: dictinary of grouped polynya keys, groups increase from 1.
- polynya_flooded: (M x N) grid of grouped polynyas, #s indicate polynya group. 0 = no polynya.

DEPENDENCIES:
import numpy as np, nump.ma as ma
import matplotlib
import matplotlib.pyplot as plt
from skimage.morphology import flood

Latest recorded update:
05-22-2024
    """
    
    # grab polynya edges 
    remain_edges = keys[polynya_edge==1]
    checked_edges = np.array([], dtype=int)

    # grid to store polynya numbers
    polynya_flooded = just_polynyas.astype(int).copy()

    # start first polynya
    polynyas = {}
    p = 1

    while len(remain_edges) > 0:

        # polynya dict
        polynyas[p] = {}
        polynyas[p]['edge'] = np.array([], dtype=int)
        polynyas[p]['inner'] = np.array([], dtype=int)

        # remove first key
        p_key = remain_edges[0]

        # record that key has been checked
        remain_edges = np.delete(remain_edges, np.where(remain_edges == p_key)[0][0], axis=None)
        checked_edges = np.append(checked_edges, p_key)

        # save polynya edge
        polynyas[p]['edge'] = np.append(polynyas[p]['edge'], p_key)

        # find where key exists
        ii, jj = np.where(keys == p_key)[0][0], np.where(keys == p_key)[1][0]

        # flood all adjacent cells (including diagonals)
        mask = flood(just_polynyas.astype(int), (ii, jj), connectivity=2)
        polynya_flooded[mask] = p

        # assign unchecked key to polynya edges or centers
        unchecked_keys = [key_ for key_ in keys[mask] if key_ not in checked_edges]

        # run through all keys and save if its on an edge or in center
        for key_ in unchecked_keys:
            if key_ in remain_edges:
                # record that key has been checked
                remain_edges = np.delete(remain_edges, np.where(remain_edges == key_)[0][0], axis=None)
                checked_edges = np.append(checked_edges, key_)
                polynyas[p]['edge'] = np.append(polynyas[p]['edge'], key_)
            else:
                polynyas[p]['inner'] = np.append(polynyas[p]['inner'], key_)

        # next polynya
        p+=1
        
        
    if show_plots != None:
        
        xx = show_plots['xx']
        yy = show_plots['yy']
        grid = show_plots['grid']
        
        # plotting
        #---------
        fig, ax = plt.subplots(figsize=(5,6))
        ax.patch.set_facecolor('lightgray')
        bounds = np.array([-0.5, 0.5, 1.5, 2.5,3.5])
        grid_norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=len(bounds)-1)
        grid_camp = matplotlib.colors.ListedColormap(['darkgray','lightsteelblue','lightgray', 'white',])
        mesh = ax.pcolormesh(xx, yy, grid, norm = grid_norm, cmap=grid_camp)
        plt.colorbar(mesh, ticks=[0,1,2,3])
        ax.pcolormesh(xx, yy, ma.masked_where(polynya_flooded==0, polynya_flooded), cmap='cmo.thermal')
        ax.pcolormesh(xx, yy, polynya_edge, cmap=matplotlib.colors.ListedColormap(['None', 'lightcoral']))


    return polynyas, polynya_flooded


def surface_type_grid(score_grid, score = 10, land_no_holes = [], show_plots = {'xx': [], 'yy': [], 'date': []}):

    """Categorize surface conditions from poylyna likelihood map.

INPUT: 
- score_grid: (M x N) grid of polynya likelihood scores
- score: score threshold to indicate polynyas (default: 10)
- land_no_holes: (M x N) boolean grid, true over land areas.
- show_plots: None if you don't want to show plots, or dict of plotting parameters 
    {'xx': (MxN) x-coordinates, 'yy': (MxN) xy-coordinates, 'date': date of data}

OUTPUT:
- grid: (M x N) category grid: 0 = land, 1 = open ocean, 2 = polynyas, 3 = ice
- water: (M x N) boolean grid of water / thin ice pixels
- just_polynyas: (M x N) boolean grid of polynya pixels
- ice: (M x N) boolean grid of ice pixels
- polynya_edge: (M x N) boolean grid of polynya edge pixels
- ocean_edge: (M x N) boolean grid of ocean edge pixels

DEPENDENCIES:
from skimage.segmentation import find_boundaries
from scipy.ndimage import binary_fill_holes
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

Latest recorded update:
05-22-2024
    """

    #=======================
    score_thresh = score
    #=======================
    
    # distinguish ice, water, from land
    water = score_grid >= score_thresh
    water = (water-land_no_holes.astype(int)) == 1
    ice = ((water == False).astype(int) - land_no_holes.astype(int)).astype(bool)
    ice_land_no_holes = binary_fill_holes(ice + land_no_holes.astype(int)).astype(int)
    just_polynyas = (ice_land_no_holes-water).astype(int)==0

    # create categorized data grid
    grid = np.zeros_like(water.astype(int))
    grid[water] = 1
    grid[just_polynyas] = 2
    grid[ice] = 3
    
    # find edges of polynya and ice
    polynya_edge = find_boundaries(just_polynyas, mode='inner')
    polynya_edges = (just_polynyas == 1) & ((find_boundaries(ice, mode='outer') - land_no_holes.astype(int))==1)
    ocean_edge = ((find_boundaries(ice, mode='outer') - land_no_holes.astype(int))==1 - polynya_edges)==1

    if show_plots != None:
        
        xx = show_plots['xx']
        yy = show_plots['yy']
        file_date = show_plots['date']
        
        # plotting
        #---------
        fig, ax = plt.subplots(figsize=(4,6))
        ax.set_title(file_date)
        ax.pcolormesh(xx, yy, score_grid, vmin=0, vmax=score_thresh, cmap='Blues')
        ax.pcolormesh(xx, yy, water, cmap=matplotlib.colors.ListedColormap(['None', 'k']))
        ax.pcolormesh(xx, yy, score_grid >= score_thresh, cmap=matplotlib.colors.ListedColormap(['None', 'indianred']))
        ax.pcolormesh(xx, yy, land_no_holes, cmap=matplotlib.colors.ListedColormap(['None', 'gray']))
        plt.show()

        # plotting
        #---------
        fig, ax = plt.subplots(figsize=(3,3))
        plt.hist(score_grid.flatten(), bins=np.arange(0,score_thresh+5,0.1))
        plt.ylim(0,2000)
        plt.vlines(score_thresh, 0,2000, colors='k', label=f'thresh = {score_thresh}')
        plt.xlabel('score')
        plt.legend()
        plt.show()

        # plotting
        #---------
        fig, ax = plt.subplots(figsize=(4,5))
        bounds = np.array([-0.5, 0.5, 1.5, 2.5,3.5])
        grid_norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=len(bounds)-1)
        grid_camp = matplotlib.colors.ListedColormap(['darkgray','lightsteelblue','dodgerblue', 'white',])
        mesh = ax.pcolormesh(xx, yy, grid, norm = grid_norm, cmap=grid_camp)
        plt.colorbar(mesh, ticks=[0,1,2,3])
        ax.pcolormesh(xx, yy, polynya_edge, cmap=matplotlib.colors.ListedColormap(['None', 'blue']))
        ax.pcolormesh(xx, yy, ocean_edge, cmap=matplotlib.colors.ListedColormap(['None', 'lightcoral']))


    return grid, water, just_polynyas, ice, polynya_edge, ocean_edge




def grab_footprint(foot, grid, ii, jj):
            
    """Grab skimage disk footprint of values around point (ii,jj) in 2D gridded data. 

INPUT: 
- foot: skimage disk footprint
- grid: (M x N) gridded data field (default: False)
- ii: (M) coordinate of data field
- jj: (N) coordinate of data field


OUTPUT:
- local_vals: list of data values within footprint of ii, jj

DEPENDENCIES:

Latest recorded update:
05-22-2024
    """
    
    # find size of footprint in either direction
    ii_size = int((foot.shape[0]-1)/2)
    jj_size = int((foot.shape[1]-1)/2)
    
    # find coordinates around edge of footprint 
    # (to handle footprints that overlap grid edges)
    ii_l = ii-ii_size
    ii_r = ii+ii_size+1
    jj_l = jj-jj_size
    jj_r = jj+jj_size+1

    # adjust bounds for footprint edge conditions
    if ii_l < 0:
        ii_l_f = -ii_l
        ii_l = 0
    else:
        ii_l_f = 0

    if jj_l < 0:
        jj_l_f = -jj_l
        jj_l = 0
    else:
        jj_l_f = 0   

    if ii_r >= grid.shape[0]:
        ii_r_f = grid.shape[0]-ii_r
        ii_r = grid.shape[0]
    else:
        ii_r_f = 1+2*ii_size

    if jj_r >= grid.shape[1]:
        jj_r_f = grid.shape[1]-jj_r
        jj_r = grid.shape[1]
    else:
        jj_r_f = 1+2*jj_size

    # crop footprint
    foot = foot[ii_l_f:ii_r_f, jj_l_f:jj_r_f]

    # crop to local values, and return those within footprint
    local_vals = grid[ii_l:ii_r, jj_l:jj_r][foot]

    return local_vals