#/////////////////////
#  download_laads ///
#///////////////////
#---------------------------------------------------------------------
# script to download NASA or VIIRS imagery given filename
#---------------------------------------------------------------------
# DEPENDENCIES:
import urllib.request
import os
#---------------------------------------------------------------------

def download_laads(file = 'MOD03.A2021060.0645.061.2021060125526.hdf', 
                   target_dir = '/Users/mackenziejewell/Desktop/temp/',
                   token_txt = '/Users/mackenziejewell/Desktop/LAADS_token.txt', 
                   allow_overwrites = False, quiet = True):

    """Use urllib to download MODIS or VIIRS imagery 

    INPUT: 
    - file: name of MODIS or VIIRS imagery file to download
        (default: 'MOD03.A2021060.0645.061.2021060125526.hdf')
    - target_dir: path to save files locally
        (default: '/Users/mackenziejewell/Desktop/temp/')
    - token_txt: path to token txt file containing NASA LAADS DAAC bearer token
        (default: '/Users/mackenziejewell/Desktop/LAADS_token.txt')
    - allow_overwrites: bool, whether or not to overwrite file if already exists locally (default: False)
    - quiet: bool, whether or not to hide print statements (default: True)

    OUTPUT: None

    DEPENDENCIES:
    import urllib.request
    import os

    Latest recorded update:
    02-07-2024
    """

    if not quiet: 
        print(f'file: {file}')
        
    # if file already exists and overwrites not allowed, don't download file
    allow_download = True
    if not allow_overwrites:
        if os.path.isfile(target_dir+file):
            allow_download = False
            if not quiet: 
                print(f'file already exists in target_dir. Will not overwrite.')

    if allow_download:
        
        # acquire bearer token from local txt
        if not quiet: 
            print(f'grab token from: {token_txt}')
        tokenfile = open(token_txt,'r')
        LAADStoken = tokenfile.read()
        tokenfile.close()

        # grab info from file name to LAADS DAAC url
        group = file.split('.')[0]
        collection = file.split('.')[3]
        # collection numbers are slightly different than those in their file names
        if collection == '061':
            collection = '61'
        elif collection == '021':
            collection = '5201'
        elif collection == '002':
            collection = '5200'
        date_string = file.split('.A')[1].split('.')[0]
        year = date_string[:4]
        jday = date_string[4:]

        # generate LAADS DAAC url
        ladsweb_url = f'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/{collection}/{group}/{year}/{jday}/{file}'
        if not quiet: 
            print(f'request >>> {ladsweb_url}')

        # request download
        opener = urllib.request.build_opener()
        opener.addheaders = [('Authorization', f'Bearer {LAADStoken}')]
        urllib.request.install_opener(opener)

        # download to target_dir
        save_file = target_dir+file
        urllib.request.urlretrieve(ladsweb_url, save_file)

        if not quiet: 
            print(f'save file to >>> {target_dir}')


