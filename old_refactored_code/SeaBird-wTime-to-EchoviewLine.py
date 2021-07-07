# -*- coding: utf-8 -*-
######################################################################
# Convert SeaBird CTD profile data to a time/depth line that can be 
# input to Echoview
# The data file has the time stamp for each measurement
# The measurements are taken at 24 Hz
# The data are not corrected - i.e., calibration coefficients are not applied
#
# The data file can be a mix of columns, with some different variables and different
# number of columns.
# The data file header string for ctd 1 & 2 in HB201805 is:
# Scan mm/dd/yyyy hh:mm:ss Latitude Longitude PrDM DepSM T090C Par Spar Sal00 Sbeox0Mg/L Sigma-e00 SvCM Flag
# Example data line:
# 1 08/13/2018 03:59:59 40.36348 -71.01316 2.694 2.694 25.6175 0.80782 1.5382 32.67 6.6833 98.43 21.3943 1533.43 0
# where columns are tab and space delimited
# The data file header string for HB201805 (#5 and up) and HB201907 is:
# Scan mm/dd/yyyy hh:mm:ss Latitude Longitude PrDM DepSM T090C Par Sal00 Sigma-e00 Sbeox0Mg/L SvCM Flag
#
# You need to delete the whacky character in Sigma-e00 (delete the "e")
#
# Export to .evl file for import to Echoview
#
# jech
#
######################################################################
import re
from datetime import datetime, timedelta
from statistics import mean

# create the path to the data
#inpath = '/net/acoustics/HB_Data/HB201805/CTD_Data'
#inpath = '/net/acoustics/HB_Data/HB201805/CTD_Data/CTD_wTime'
#!! change the scan line (~line 75) for 2018 or 2019 data !!#
#inpath = '/net/acoustics/HB_Data/HB201805/CTD_Data/CTD_wTime'
#inpath = '/net/acoustics/HB_Data/HB201907/CTD_Data'
#inpath = '/home/jjech/Desktop/CTD_Data'
#inpath = '/media/mjech/Mac Passport/GU201905/CTD_Data'
inpath = '/Volumes/GeringSSD/GU201905_CTD'

# toffset is the time offset in seconds for the cast relative to 
# the acoustic data
toffset = 0   # seconds

# the file names are in a file (bfile) with their respective time offsets
bfile = inpath+'/'+'CTDtoEVL.list'
# read the entire file as a list of lines
with open(bfile, 'r') as infile:
    CTD_files = infile.readlines()

# cycle through the files
for j in range(0, len(CTD_files)):
    # the offset in time between the FS70 data and real/actual time in seconds
    (infn, toffset) = CTD_files[j].split(',')
    infn = inpath+'/'+infn.strip()
    toffset = int(toffset)
    print('Doing: ', infn, 'with toffset: ', toffset)

    # read the entire file as a list of lines
    # added extra encoding parameter to avoid problem with reading whacky character in Sigma-e00 
    # https://stackoverflow.com/questions/19699367/for-line-in-results-in-unicodedecodeerror-utf-8-codec-cant-decode-byte
    with open(infn, 'r', encoding = "ISO-8859-1") as infile:
        CTD_data = infile.readlines()

    # create the output file using the inputfile name
    outfn = infn.replace('.asc', '.evl')
    outfile = open(outfn, 'w')
    # the first output header line is EVBD 3, then the Echoview version
    headerline = 'EVBD 3 9.0.298.34146'
    outfile.write(headerline+'\n')

    # read the data file in. create a dictionary with date and time as the key 
    # and depth as the value. The data have time to the second, so all depths
    # within a second are kept as a list, then take the mean depth for each
    # second
    data = {}
    for line in CTD_data[1:]:
        # strip the end of line
        line.rstrip()
        
        #(Scan, date, time, lat, lon, PrDM, depth, ToC, Par, Spar, sal, SbeoxMG, \
        #SbeoxPS, sigma, Sv, flag) = line.split()
        #(Scan, date, time, lat, lon, PrDM, depth, ToC, Par, sal, SbeoxMG, \
        # SbeoxPS, sigma, Sv, flag) = line.split()
        # the variables change among files so just split and fortunately the date, time,
        # lat, lon, PrDM, and depth seem to be consistent among all files
        (Scan, date, time, lat, lon, PrDM, depth, tmpvars) = line.split(maxsplit=7)
#        date = tmp[1]
#        time = tmp[2]
#        depth = tmp[6]
        # create a dictionary based on date & time
        tmpkey = date+'-'+time
        if tmpkey in data:
            data[tmpkey].append(float(depth))
        else:
            data[tmpkey] = [float(depth)]
    
    # the number of lines will be the number of keys
    outfile.write(str(len(data))+'\n')
    # need the dict in chronological order
    # take the mean of the depths for each second
    for k, v in sorted(data.items()):
        depth_mean = mean(v)
        (date, time) = str.split(k, sep='-')
        (month, day, year) = str.split(date, sep='/')
        (hour, minute, second) = str.split(time, sep=':')
        #print(year+month+day,' ', hour+minute+sec+'0000 ', round(depth_mean, 1), ' 3\n')
        # need to accomodate a time offset 
        tmptime = datetime(int(year), int(month), int(day), \
                  int(hour), int(minute), int(second))+timedelta(seconds=toffset)
        outfile.write(tmptime.strftime('%Y%m%d %H%M%S'+'0000 ')+\
                      str(round(depth_mean, 1))+' 3\n')
                      
    # clean house
    outfile.close()
### end of main

# the file name
#infn = inpath+'/'+'Bigelow1805_CTD_Cast8.asc'
#toffset = 5*60+53  # CTD 8 HB201805
#infn = inpath+'/'+'CTD008.asc'
#toffset = 3*60-3  # CTD 8 HB201805

#infn = inpath+'/'+'Bigelow1805_CTD_Cast5.asc'
#toffset = 5*60+50  # CTD 5 HB201805
#infn = inpath+'/'+'CTD005.asc'
#toffset = 3*60  # CTD 5 HB201805

#infn = inpath+'/'+'Bigelow1805_CTD_Cast6.asc'
#toffset = 5*60+50  # CTD 6 HB201805
#infn = inpath+'/'+'CTD006.asc'
#toffset = 3*60+35  # CTD 6 HB201805

#infn = inpath+'/'+'Bigelow1805_CTD_Cast7.asc'
#toffset = 5*60+50  # CTD 7 HB201805
#infn = inpath+'/'+'CTD007.asc'
#toffset = 3*60  # CTD 7 HB201805

#infn = inpath+'/'+'Bigelow1805_CTD_Cast9.asc'
#toffset = 4*60+35  # CTD 9 HB201805
#infn = inpath+'/'+'CTD009.asc'
#toffset = 2*60+10 # CTD 9 HB201805

#infn = inpath+'/'+'Bigelow1805_CTD_Cast10.asc'
#toffset = 5*60+5  # CTD 10 HB201805
#infn = inpath+'/'+'CTD010.asc'
#toffset = 3*60+13  # CTD 10 HB201805

#infn = inpath+'/'+'ctd001.asc'
#toffset = 3*60+50  # CTD 1 HB201805

#infn = inpath+'/'+'ctd002.asc'
#toffset = 3*60  # CTD 2 HB201805

