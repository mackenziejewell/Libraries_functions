#/////////////////////
#  write_nc_file  ///
#///////////////////
#---------------------------------------------------------------------
# Write xarray dataset to NETCDF4 file 
# after checking whether file already exists.
#---------------------------------------------------------------------
#//////////////////////
#  write_csv_file  ///
#////////////////////
#---------------------------------------------------------------------
# Write xarray dataset to csv file 
# after checking whether file already exists.
#---------------------------------------------------------------------









#/////////////////////
#  write_nc_file  ///
#///////////////////
#---------------------------------------------------------------------
# Write xarray dataset to NETCDF4 file 
# after checking whether file already exists.
#---------------------------------------------------------------------
# DEPENDENCIES:
import os
import xarray
#---------------------------------------------------------------------
    
def write_nc_file(ds, outfile_path = 'outfile.nc', allow_overwrites = False):
    
    """Write xarray dataset to NETCDF4 file after checking whether file already exists.

INPUT: 
- ds: xarray dataset to be saved
- outfile_path: location to save, ending in '.nc' (default: 'outfile.nc' in local folder)
- allow_overwrites: bool, whether or not to overwrite if file already exists 
  (default: False)

DEPENDENCIES:
import os
import xarray

Latest recorded update:
08-18-2022
    """

    # check whether file is being saved with nc suffix, as given by outfile_path
    assert str(outfile_path[-3:])=='.nc', f"outfile should end in '.nc', not '{outfile_path[-3:]}'"
    # check whether file was given as xarray dataset
    assert str(type(ds)) == "<class 'xarray.core.dataset.Dataset'>", f"ds should be xarray dataset, not {str(type(ds))}"
    
        
    # check whether the txt file already exists in the Destination
    #-------------------------------------------------------------
    isFile = os.path.isfile(outfile_path)
    
    # if file already exists, check if it should be overwritten
    #----------------------------------------------------------
    if isFile == True:
        print('\n!!! file already exists !!!\n')
        print(f'>>> {outfile_path}')
        # determine whether or not to write file
        if allow_overwrites == True:
            print('\n>>> Overwriting existing file.')
            WriteFile = True
        else:
            print('\n>>> Will not overwrite.')
            WriteFile = False
    else:
        WriteFile = True

    # write file if allowed
    #----------------------
    if WriteFile == True:
        print(f'Saving {outfile_path}')
        ds.to_netcdf(path=outfile_path, mode='w', format="NETCDF4")
        
        
        

#//////////////////////
#  write_csv_file  ///
#////////////////////
#---------------------------------------------------------------------
# Write xarray dataset to csv file 
# after checking whether file already exists.
#---------------------------------------------------------------------
# DEPENDENCIES:
import os
import xarray
#---------------------------------------------------------------------
def write_csv_file(ds, outfile_path = 'outfile.csv', allow_overwrites = False):
    
    """Write xarray dataset to NETCDF4 file after checking whether file already exists.

INPUT: 
- ds: xarray dataset to be saved
- outfile_path: location to save, ending in '.csv' (default: 'outfile.csv' in local folder)
- allow_overwrites: bool, whether or not to overwrite if file already exists 
  (default: False)

DEPENDENCIES:
import os
import xarray

Latest recorded update:
04-20-2022
    """
    
    # check whether the txt file already exists in the Destination
    #-------------------------------------------------------------
    isFile = os.path.isfile(outfile_path)
    # if file already exists, check if it should be overwritten
    #----------------------------------------------------------
    if isFile == True:
        print('\n!!! file already exists !!!\n')
        print(f'>>> {outfile_path}')
        # determine whether or not to write file
        if allow_overwrites == True:
            print('\n>>> Overwriting existing file.')
            WriteFile = True
        else:
            print('\n>>> Will not overwrite.')
            WriteFile = False
    else:
        WriteFile = True

    # write file if allowed
    #----------------------
    if WriteFile == True:
        print(f'Saving {outfile_path}')
        ds.to_csv(outfile_path, index=False)