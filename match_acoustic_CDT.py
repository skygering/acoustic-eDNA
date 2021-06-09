from plot_CTD_data import evl_list
from datetime import datetime, timedelta
import glob
import os

path = '/Volumes/GeringSSD/GU201905_CTD' # read in data - CHANGE TO MATCH LOCAL FILE STRUCTURE
# Use CTDtoEVL list to plot - assumes that you just ran SeaBird-wTime-toEchoviewLine.py
ctd_list_file = path+'/'+'CTDtoEVL.list'

evl_files = evl_list(ctd_list_file) # get list of evl files

evl_dates_dic={}
for file in evl_files:
    with open(file, 'r', encoding = "ISO-8859-1") as infile:
        evl_data = infile.readlines()
    start = evl_data[2] # first line of data after header
    start.rsplit()
    (sd, st, _, _) = start.split(maxsplit=3)
    start_dt = datetime.strptime(sd+st, "%Y%m%d%H%M%S%f") # start time of CDT data
    end = evl_data[-1]
    end.rsplit()
    (ed, et, _, _) = start.split(maxsplit=3)
    end_dt = datetime.strptime(ed+et, "%Y%m%d%H%M%S%f") # end time of CTD data
    evl_dates_dic[os.path.basename(file)] = (start_dt, end_dt) # dictionary of files with tuple start and end time

raw_list = glob.glob("/Volumes/GeringSSD/GU1905_Acoustic/EK60/*.raw") # Read in all .raw files from directory

date_raw_dic={}
for filepath in raw_list:
    file = os.path.basename(filepath)
    cruise, date, time = file.replace(".raw", "").split("-") # Break apart .raw filename to get cruise, date, and time
    start_dt = datetime.strptime(date[1:] + time[1:], "%Y%m%d%H%M%S") # Start of acoustic data collection for file
    date_raw_dic[start_dt] = file # Dictionary of start times and .raw files

outfile = open(path + "/evl_raw_matches.list", 'w')
headerline = 'Cruise ' + cruise + " CTD and matching acoustic (.raw) files"
outfile.write(headerline+'\n')

start_delay = timedelta(minutes=15) # Time delay for starting CTD data collection after acoustics 
for evl_file in evl_dates_dic:
    (stime, etime) = evl_dates_dic[evl_file]
    raw_list = []
    for date in date_raw_dic:
        # If acoustic start is between CDT times or start_delay before, record that raw file
        if (date > stime and date < etime) or abs(date - stime) < start_delay:
            raw_list.append(date_raw_dic[date])
    outfile.write(evl_file + " ".join(raw_list) + "\n") # write CDT file (.evl) followed by overlapping .raw files

outfile.close()


