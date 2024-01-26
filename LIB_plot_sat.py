#/////////////////////////////
#  check_imagery_matches  ///
#///////////////////////////
#---------------------------------------------------------------------
# Confirm matches for downloaded geo and imagery files
#---------------------------------------------------------------------
#////////////////////////
#  pair_images_meta  ///
#//////////////////////
#---------------------------------------------------------------------
# Load metadata of MODIS image files from folder.
#---------------------------------------------------------------------
#///////////////////////////
#  plot_singleband_sat  ///
#/////////////////////////
#---------------------------------------------------------------------
# Make plot of projected single-band MODIS/VIIRS/etc image from imagery and geo file.
#---------------------------------------------------------------------


#/////////////////////////////
#  check_imagery_matches  ///
#///////////////////////////
#---------------------------------------------------------------------
# Confirm matches for downloaded geo and imagery files
#---------------------------------------------------------------------
# DEPENDENCIES
#---------------------------------------------------------------------
def check_imagery_matches(geo_file_list = [], img_file_list = [], quiet = True):
    
    """Confirm that each geo file has a corresponding imagery match and vice versa. Works for folders containing MODIS/terra, MODIS/aqua, VIIRS/JPSS1, VIIRS/NPP files downloaded frmo NASA LAADS DAAC.
    
INPUT:
- geo_file_list: list of geolocation files, beginning in 
'VJ103MOD', 'VNP03MOD', 'MOD03', or 'MYD03' depending on the platform.
- img_file_list: list of imagery files, beginning in 
'VJ102MOD', 'VNP02MOD', 'MOD021KM', or 'MYD021KM' depending on the platform.
- quiet: bool, whether or not to suppress prints (default: True)

OUTPUT:
- missing_geo: list of missing geo files, given the downloaded img_file_list
- missing_img: list of missing imagery files, given the downloaded geo_file_list

DEPENDENCIES:

Latest recorded update:
07-20-2023
    """
    
    # lists to store missing geo and img files
    missing_geo = []
    missing_img = []
    
    
    # check for missing imagery files
    #--------------------------------
    for geo_file in geo_file_list:
    
        # grab platform info and date info
        platform = geo_file[:geo_file.find('A')-1]
        date_info = geo_file[geo_file.find('A'):geo_file.find('A')+17]
        
        # find corresponding platform info
        if str(platform) == 'VJ103MOD':
            corr_plat = 'VJ102MOD'
        elif str(platform) == 'VNP03MOD':
            corr_plat = 'VNP02MOD'
        elif str(platform) == 'MOD03':
            corr_plat = 'MOD021KM'
        elif str(platform) == 'MYD03':
            corr_plat = 'MYD021KM'
        else:
            print(f' >> ! unrecognized platform {platform} in geo_file {geo_file}')
        
        # construct matching image file string portion
        file_match = corr_plat+'.'+date_info

        # look for img matching geo_file in list, 
        # else print that file is missing and save
        # to missing list
        match = False
        for img_file in img_file_list:
            if file_match in img_file:
                match = True
                break
        if not match:
            missing_img.append(file_match)
            if not quiet:
                print(f'- geo_file {geo_file} missing image file')

    # print missing imgs if any are missing
    if len(missing_img) > 0:
        if not quiet:
            print(f'\n{len(missing_img)} missing imagery files\n-------------------------')
            for img in missing_img:
                print(img)
        print()
    
    # check for missing imagery files
    #--------------------------------
    for img_file in img_file_list:
    
        # grab platform info and date info
        platform = img_file[:img_file.find('A')-1]
        date_info = img_file[img_file.find('A'):img_file.find('A')+17]
        
        # find corresponding platform info
        if str(platform) == 'VJ102MOD':
            corr_plat = 'VJ103MOD'
        elif str(platform) == 'VNP02MOD':
            corr_plat = 'VNP03MOD'
        elif str(platform) == 'MOD021KM':
            corr_plat = 'MOD03'
        elif str(platform) == 'MYD021KM':
            corr_plat = 'MYD03'
        else:
            print(f' >> ! unrecognized platform {platform} in img_file {img_file}')
        
        # construct matching image file string portion
        file_match = corr_plat+'.'+date_info

        # look for img matching geo_file in list, else print 
        # that file is missing and save to missing list
        match = False
        for geo_file in geo_file_list:
            if file_match in geo_file:
                match = True
                break
        if not match:
            missing_geo.append(file_match)
            if not quiet:
                print(f'- img_file {img_file} missing geo file')

    # print missing geo if any are missing
    if len(missing_geo) > 0:
        if not quiet:
            print(f'\n{len(missing_geo)} missing geo files\n-------------------------')
            for geo in missing_geo:
                print(geo)
        print()
    
    return missing_geo, missing_img




#////////////////////////
#  pair_images_meta  ///
#//////////////////////
#---------------------------------------------------------------------
# Load metadata of MODIS image files from folder.
#---------------------------------------------------------------------
# DEPENDENCIES
import glob
import os
import numpy as np
# homemade: 
from LIB_plot_VIIRS import get_VIIRS_date
from LIB_plot_MODIS import get_MODISdate
from LIB_plot_sat import check_imagery_matches
#---------------------------------------------------------------------

def pair_images_meta(MainFolder = [], SingleFolder = [], 
                     sensor = 'MODIS', 
                     satellite_labels = [('MOD03','MOD021KM'), ('MYD03','MYD021KM')], 
                     min_geofile_sizeMB = [], min_imfile_sizeMB = [],
                     max_diff_minutes = 20):
    
    """Load metadata of MODIS image files from folder. Pair geo and image and group by date.
    Does not open files, all operations done from file names so make sure
    the MODIS files are saved with the original names as downloaded from the LAADS DAAC
    
INPUT:
- MainFolder:main directory where subdirectories (one step down) contains images to check(default: []) 
- SingleFolder: directory where images are stored (default: []) 
- sensor: string describing sensor from satellite data (currently either 'VIIRS' for suomiNPP (.nc files) or 'MODIS' (.hdf) from NASA LAADS DAAC
- satellite_labels: list of tuples for geolocation and imagery tags for each satellite surce
                    (default: [('MOD03','MOD021KM'), ('MYD03','MYD021KM')] for Terra/MODIS, Aqua/MODIS)
- min_geofile_sizeMB: minimum accepted geofile size (in MB) without raising error
                      used to check for corrupted files which have too small of file size (normally > 28 MB)
                      (default: [] in which case it won't check file size)
- min_imfile_sizeMB: minimum accepted image file size (in MB) without raising error
                     used to check for corrupted files which have too small of file size (normally > 55 MB)
                     (default: [] in which case it won't check file size)
- max_diff_minutes: maximum time difference between images to use for pairing
                    (default: 20)

OUTPUT:
- Image_Meta_paired: M x 5 array of paired image metadata 
                     [date, geo_filename, image_filename, filepath, pair_index]
                     also paired images/dates are printed

DEPENDENCIES:
import glob
import os
import numpy as np
# homemade: 
from LIB_plot_VIIRS import get_VIIRS_date
from LIB_plot_MODIS import get_MODISdate
from LIB_plot_sat import check_imagery_matches

Latest recorded update:
07-20-2023

    """

    # assertions for input data
    #--------------------------
    assert sensor in ['MODIS', 'VIIRS'], f"Unrecognized satellite type, got: {sensor}"
    if str(sensor) == 'VIIRS':
        file_type = ".nc"
    elif str(sensor) == 'MODIS':
        file_type = ".hdf"
        
    Image_Meta = []
    
    # find all folders in main directory
    # or look in single folder provided
    if MainFolder != []:
        assert os.path.exists(MainFolder), f"{MainFolder} not an existing directory"
        folder_list = glob.glob(MainFolder+"*/")
        print('Search within main folder: {}'.format(MainFolder))
    elif SingleFolder != []:
        assert os.path.exists(SingleFolder), f"{SingleFolder} not an existing directory"
        folder_list = glob.glob(SingleFolder)
        print('Search in single folder: {}'.format(SingleFolder)) 
                  
    # find list of hdf/nc files in each folder
    #--------------------------------------
    for folder in folder_list:
        
        # grab all hdf/nc files in folder and start empty
        # lists to fill with geo and image files
        #---------------------------------------------
        file_list = sorted(list(glob.glob1(folder, f"*{file_type}")));

        image_list,geo_list = [],[]
        # check for an even number of hdf/nc or nc files in folder
        # since num geo and image files should match
        #-----------------------------------------------
        if len(file_list)%2!=0:
            print(f'Odd number of {file_type} files found in folder. Should be one geo file per imagery file.')
#             break
        # generate lists of geolocation and imagery files
        # run through filenames in file_list
        # for each type of satellite in satellite_labels:
        # if geolocation label in filename, add to geo_list
        # if imagery label in filename, add to geo_list
        #------------------------------------------------
        for file in file_list:  
            # grab size (MB) of file
            if min_geofile_sizeMB != [] and min_imfile_sizeMB != []:
                file_sizeMB = os.path.getsize(folder+file)/(1000**2) 
            for satellite in satellite_labels:
                if satellite[0] in file:
                    if min_geofile_sizeMB != []:
                        # if file size is less than given min_geofile_sizeMB, throw error since file
                        # was probably corrupted upon download. Usually 30-40 MB
                        if file_sizeMB < min_geofile_sizeMB:
                            print('{}\n{}\n{}'.format('='*len('POSSIBLE ERROR:'),'POSSIBLE ERROR:','='*len('POSSIBLE ERROR:')))
                            print('File {} is only {:.1f} MB, likely corrupted'.format(file,file_sizeMB))
                            break
                    geo_list.append(file)
                if satellite[1] in file:
                    if min_imfile_sizeMB != []:
                        # if file size is less than 55 MB, throw error since file
                        # was probably corrupted upon download. Usually 60+ MB
                        if file_sizeMB < min_imfile_sizeMB:
                            print('{}\n{}\n{}'.format('='*len('POSSIBLE ERROR:'),'POSSIBLE ERROR:','='*len('POSSIBLE ERROR:')))
                            print('File {} is only {:.1f} MB, likely corrupted\n'.format(file,file_sizeMB))
                            break
                    image_list.append(file)
                    
         
        # check for imagery and geo file matches
        #-------------------------------------------
        mg, mi = check_imagery_matches(geo_file_list=geo_list, img_file_list=image_list, quiet=False)
        
        # for each image, create list 
        # [date, geo_filename, image_filename, filepath]
        while len(geo_list)>0:
            # set current geo_file to first file in list
            # and delete geo_file from geo_list
            # and grab date from geo_file
            #-------------------------------------------
            geo_file = geo_list.pop(0)  
            
            if str(sensor) == 'VIIRS':
                ImageDate = get_VIIRS_date(geo_file)
            elif str(sensor) == 'MODIS':
                ImageDate = get_MODISdate(geo_file)
                
            # determine satellite source of geo_file
            #---------------------------------------
            for ii in range(len(satellite_labels)):
                if satellite_labels[ii][0] in geo_file:
                    satellite_image_source = satellite_labels[ii][1]
                    
            # search for matching imagery file
            # and delete it from image_list
            #---------------------------------
            ii = 0
            while len(image_list)>0: 
                # print error if search gets to end of list before finding a match
                if ii==len(image_list):
                    print(f'Corresponding image match could not be found for geo_file {geo_file}')
                    break
                # if satellite source matches geofile
                # check if imagery date matches geo_file
                # delete image_file from geo_list
                #------------------------------------
                if satellite_image_source in image_list[ii]:
                    if str(sensor) == 'VIIRS':
                        if get_VIIRS_date(image_list[ii])==ImageDate:
                            image_file = image_list.pop(ii)
                            break
                    elif str(sensor) == 'MODIS':
                        if get_MODISdate(image_list[ii])==ImageDate:
                            image_file = image_list.pop(ii)
                            break
                ii+=1
            Image_Meta.append([ImageDate, geo_file, image_file, folder])
    Image_Meta = sorted(Image_Meta)

    
    # FIND IMAGE PAIRS AND ADD PAIR_INDEX TO META
    #--------------------------------------------
    Image_Meta_paired = np.array([])
    pair_index = 0
    while len(Image_Meta)>0: 
        
        # pull current first image out of Image_Meta list
        # add a pair index to its metadata and add to Image_Meta_copy
        current_image = Image_Meta.pop(0)
        current_image.append(pair_index)
        Image_Meta_paired = np.append(Image_Meta_paired,current_image)
    
        
        # run through rest of images to search for pair
        # find difference in time between current image remaining images
        # if within max_diff_minutes mins, add to Image_Meta_copy with same pair_index
        indices_to_pair = []
        for check_image in Image_Meta:
            time_diff = np.abs(current_image[0]-check_image[0]).total_seconds()/60
            if time_diff <= max_diff_minutes: 
                indices_to_pair.append(Image_Meta.index(check_image))
        
        images_to_add_to_pair = np.array([])
        # in reverse order, pop paired values from Image_Meta
        for current_index in indices_to_pair[::-1]:
            pair_image = Image_Meta.pop(current_index)
            pair_image.append(pair_index)
            images_to_add_to_pair = np.append(images_to_add_to_pair, pair_image)
        images_to_add_to_pair = np.reshape(images_to_add_to_pair, [int(len(indices_to_pair)),5])
        
        # in original sequential order, add to to Image_Meta_paired with correct index
        for pair_image in images_to_add_to_pair[::-1]:
            Image_Meta_paired = np.append(Image_Meta_paired,pair_image)
                
        # move on to next of remaining images
        pair_index+=1
        
    Image_Meta_paired = np.reshape(Image_Meta_paired,[int(len(Image_Meta_paired)/5),5])

    # print image dates by pairs
    #---------------------------
    num_pairs = np.max(Image_Meta_paired[:,4])+1
    for ii in range(num_pairs):
        print('Pair {}\n------'.format(ii))
        for jj in np.where(Image_Meta_paired[:,4] ==ii)[0]:
            print(Image_Meta_paired[jj,0])
        print()
        
        
    return Image_Meta_paired





#///////////////////////////
#  plot_singleband_sat  ///
#/////////////////////////
#---------------------------------------------------------------------
# Make plot of projected single-band MODIS/VIIRS/etc image from imagery and geo file.
#---------------------------------------------------------------------
# DEPENDENCIES:
#-------------
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors                                
#---------------------------------------------------------------------

def plot_singleband_sat(ax, LAT, LON, _image_, SourceCCRS = ccrs.PlateCarree(), 
                        cmap = 'Greys',  cscale = [1.8, 6.8], shading = 'gouraud',zorder=1): 
    
    
    """Make plot of projected single-band satellite image from imagery and geo file.
    
INPUT:
- ax : figure axis
- LAT, LON: lat/lon grids of image to plot (or list of grids if overlaying multiple images simultaneously)
- _image_ : image to plot (or list of images if overlaying multiple images simultaneously)  
- SourceCCRS: CRS of input data, should use ccrs.PlateCarree() (default: ccrs.PlateCarree())
- cmap: colormap to use (default: 'Greys')
- cscale: normalization for data in colormap (default: [min, max] = [1.8, 6.8])
            Mx1 lists of floats (M=2,3,4) either as:
              - [vmin, vmax] 
              - [vmin, midpoint, vmax] 
              - [vmin, midpoint1, midpoint2, vmax]
- shading: pcolormesh shading style (default: 'gouraud')
- zorder: drawing order of layer (default: 1)


OUTPUT:
 adds projected single-band satellite image to plot

DEPENDENCIES:
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors

Latest recorded update:
10-07-2022

    """
    # colormap normalization
    #-----------------------
    # class for colormap normalizing scaling options
    class TwopointNormalize(matplotlib.colors.Normalize):
        def __init__(self, vmin=None, vmax=None, vmid1=None, vmid2=None, clip=False):
            self.vmid1 = vmid1
            self.vmid2 = vmid2
            super().__init__(vmin, vmax, clip)
        def __call__(self, value, clip=None):
            x, y = [self.vmin, self.vmid1, self.vmid2, self.vmax], [0, 0.33,0.66, 1]
            return np.ma.masked_array(np.interp(value, x, y))

    # determine color normalization
    if len(cscale)==2:
        norm=matplotlib.colors.Normalize(vmin=cscale[0], vmax=cscale[1])
    elif len(cscale)==3:
        norm=matplotlib.colors.TwoSlopeNorm(vmin=cscale[0], vcenter=cscale[1],vmax=cscale[2])
    else:
        norm = TwopointNormalize(vmin=cscale[0], 
                                        vmid1=cscale[1], 
                                        vmid2=cscale[2], 
                                        vmax=cscale[3])
        
    # plot level1b image
    #-------------------
    
    # for each file in _image_ list
    for ii in range(0,len(_image_)):
        ax.pcolormesh(LON[ii], LAT[ii], _image_[ii], cmap = cmap, norm = norm,
                      transform=SourceCCRS, zorder = zorder, shading = shading)
   # use 'gouraud' shading in pcolormesh since it doesn't like that it's given cell centers rather than corners
        
    # remove automatic image border
    ax.outline_patch.set_linewidth(0)