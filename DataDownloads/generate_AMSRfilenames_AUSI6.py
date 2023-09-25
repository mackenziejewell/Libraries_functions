import os
import calendar

# list all years of data to download
#===================================
all_years = [2013]
#===================================

# list all month numbers to download each year
#=============================================
all_months = [1,2,3,4,5,6,10,11,12]
#=============================================

# name of text file to create
#====================================
out_txt_name = "2_AMSRfilenames.txt"
#====================================

# structure of https filenames / directory structures
#----------------------------------------------------
# {} contain date information, reference:
# example file name for data from Jan 1, 2013
# 'https://n5eil01u.ecs.nsidc.org/AMSA/AU_SI6.001/2013.01.01/AMSR_U2_L3_SeaIce6km_B04_20130101.he5'
string_base = 'https://n5eil01u.ecs.nsidc.org/AMSA/AU_SI6.001/{}.{}.{}/AMSR_U2_L3_SeaIce6km_B04_{}{}{}.he5'

all_files = []

# all years 
#----------
for year in all_years:

    # all months 
    #-----------
    for month in all_months:

        # find number of days in given month
        (m, days) = calendar.monthrange(year, month)

        # all days 
        #---------
        # run through each day int
        for dd in range(days):

            # generate day int, starting from day = 1
            day = dd + 1

            # convert to 0-led strings
            mon_s = str(month).zfill(2)
            day_s = str(day).zfill(2)

            # create filename and save to list
            file = string_base.format(year, mon_s, day_s, year, mon_s, day_s)
            all_files.append(file)

# save generated filenames to txt file
#-------------------------------------
with open(out_txt_name, "w") as a:
    for file in all_files:
        a.write(str(file) + os.linesep)


print(f'script to generate filenames for NSIDC AU_SI6 (doi: 10.5067/NX1R09ORNOZN)\n')
print(f'>>> years: {all_years}')
print(f'>>> months: {all_months}')
print(f'>>> create txt file: "{out_txt_name}" containing {len(all_files)} file names')