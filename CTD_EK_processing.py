from echolab2.processing import line
from datetime import date, datetime, timedelta
from statistics import mean
import numpy as np
import os
import math
import glob

def asc_to_evl(ctd_list_fn, infile_path, outfile_path, echoview_version = "EVBD 3 9.0.298.34146"):
    '''
    asc_to_evl: takes a file with a list of .asc files and time offsets for cast relative to acoustic data and writes .evl files;
                each .evl file has sample times, mean depth, and status of the data (3 if good, 2 if depth is negative of NaN)
   
    Inputs: 
            ctd_list_fn (string) - file name (suggested name 'CTDtoEVL.list') for a file which has the following format with example below:

                Format: example_ctd.asc, toffset (** Note this is for one line of the file **)
                Example: ctd001.asc, 0

                where example_ctd.asc is the name of the CTD cast and toffset is the time offset in seconds for that cast 
                relative to the acoustic data. Repeat the format on each line of the file for each cast

            infile_path (string) - path to all .asc files and path to 'CTDtoEVL.list'
            outfile_path (string) - path where you want all .evl files saved to
            echoview_verion (string) - optional input denoting the version of echoview you're using if you want to read samples into Echoview
                
                Format: EVBD 3 echoview_version
                Example: EVBD 3 9.0.298.34146

    Outputs: a .evl file for each .asc file saved at outfile_path that takes into account the time offset of cast and acoustics
             and denotes good/bad data points with 3 or 2 respectivly. It has the following format, with one line per data point:

                Format: sample_date sample_time sample_mean_depth sample_status
                Example: 20191018 1142460000 3.0 3
            
            Note there are also 2 headers, needed for reading into Echoview: an Echoview version and the number of data points
    '''
    
    ctd_list_path = os.path.normpath(infile_path + '/'+ ctd_list_fn) # list is assumed to exist

    # read the entire file as a list of lines
    with open(ctd_list_path, 'r') as infile:
        ctd_list = infile.readlines()
    
    # cycle through file/time offset list
    for line in ctd_list:
        (asc_infn, toffset) = line.split(',') # each line has .asc file and time offset between cast and acoustic data
        asc_infn_path = os.path.normpath(infile_path + '/' + asc_infn) # location of all .asc files
        toffset = int(toffset)

        print('Doing: ', asc_infn_path, 'with toffset: ', toffset)

        # open each .asc file and read data
        with open(asc_infn_path, 'r', encoding = "ISO-8859-1") as asc_file: # encoding allows reading of epsilons in header
            ctd_data = asc_file.readlines()

        data = {}
        for sample in ctd_data[1:]: # first line is header for each column
            sample.rstrip()
            # there are more than 7 data types, but we only need data/time/depth
            (Scan, date, time, lat, lon, PrDM, depth, tmpvars) = sample.split(maxsplit=7)
            (month, day, year) = str.split(date, sep='/')
            (hour, minute, second) = str.split(time, sep=':')
            # create date time object for each sample
            sample_dt = datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))+timedelta(seconds=toffset)

            # same the depths recorded for each time in a dictionary
            if sample_dt in data:
                data[sample_dt].append(float(depth)) # more than one depth saved for a single time
            else:
                data[sample_dt] = [float(depth)]

        # Preparing outfile for writing
        evl_outfn = os.path.normpath(outfile_path + "/" + asc_infn.strip().replace("asc", "evl"))
        outfile = open(evl_outfn, 'w')
        headerline = echoview_version # the first output header line is EVBD 3, then the Echoview version
        outfile.write(headerline+'\n')
        outfile.write(str(len(data))+'\n') # second header is the number of data points - number of time stamps

        for sample_dt, sample_depths in sorted(data.items()):
            depth_mean = mean(sample_depths) # mean depth for any one time
            sample_status = 3 # good data
            if depth_mean < 0 or depth_mean is math.nan:
                sample_status = 2 # bad data
            outfile.write(sample_dt.strftime('%Y%m%d %H%M%S'+'0000 ')+\
                      str(round(depth_mean, 1))+ ' ' + str(sample_status) + '\n')
        
        outfile.close() # close outfile
            

def asc_from_list(ctd_list_fn):
    '''
    asc_from_list: takes a file with a list of .asc files and time offsets for cast relative to acoustic data and
                   returns a list of all of the .asc files
    Inputs: ctd_list_fn (string) - file name (suggested name 'CTDtoEVL.list') for a file which has the following format with example below:

                Format: example_ctd.asc, toffset (** Note this is for one line of the file **)
                Example: ctd001.asc, 0

                where example_ctd.asc is the name of the CTD cast and toffset is the time offset in seconds for that cast 
                relative to the acoustic data. Repeat the format on each line of the file for each cast.
            *This is a file needed for asc to evl - this function is to make further use of that file
    Outputs: a string list of .asc file names from ctd_list_fn
    '''
    ctd_list_fn = os.path.normpath(ctd_list_fn)
    with open(ctd_list_fn, 'r') as infile:
        ctd_files = infile.readlines() # read in each CDT file name and the toffset
    asc_files = []
    for ctd in ctd_files:
        # Get the filename and replace .asc with .evl - you now have all .evl filenames
        (file, _) = ctd.split(',')
        asc_files.append(file)
    return asc_files

def cast_new_extension(ctd_file_list, old_extension, new_extension, descriptor=""):
    """
    cast_new_extension: takes a list of filenames with the same extensions and replaces the extensions with new extensions;
                        can also take an optional file descriptor to add to the filename before the new extension
    Inputs: ctd_file_list (string list) - list of sting filenames with all the same extensions
            old_extension (string) - extension of files in ctd_file_list
            new_extension (string) - new extension for all filenames in ctd_file_list
            descriptor (string) - optional input
    Outputs: (string list) - list of string filenames with new extension and optional descriptor
            Example: ctd001.evl becomes ctd001_descriptor.jpg where old_extension = ".evl" and new_extension = ".jpg"
    """

    new_files = []
    for file in ctd_file_list:
        file = file.strip().replace(old_extension, "")
        if len(descriptor) != 0:
            file = file + "_" + descriptor
        file += new_extension
        new_files.append(file)
    return new_files

def match_raw_evl(evl_files_list, raw_files_list, outfile_path = "", evl_inpath = "", raw_inpath = ""):
    '''
    match_raw_evl: matches acoustic (raw) files to CTD casts (evl)file by time
    Inputs: evl_files_list (string list) - list of .evl filenames - if they don't have a path within the filename
            then the evl_inpath variable is needed
            raw_files_list (string list) - list of .raw filenames - if they don't have a path within the filename
            then the raw_inpath variable is needed
            outfile_path (string) - path where the list of pairs of evl files and raw files is printed
            evl_inpath (string) - path to evl files, ONLY use if the filenames in evl_files_list don't have a path
            raw_inpath Istring) - path to raw files, ONLY use if the filenames in raw_files_list don't have a path
    Outputs: evl_raw_dic (string dictionary) - keys are .evl files and the values are lists of .raw files that overlap
             with that .evl file within a 30 minute time delay

             if outfile_path is not empty - file with list of .evl files and .raw files 

                Format: example_cast.evl raw_overlap_1.raw raw_overlap_2.raw
                Example: ctd017.evl GU19_05-D20191027-T140551.raw GU19_05-D20191027-T143423.raw

             There is one line per CTD cast
    '''
    evl_dates_dic={}
    for file in evl_files_list:
        evl_line = line.read_evl(os.path.normpath(evl_inpath+ "/" + file))
        date_time = evl_line.ping_time
        date_time = np.sort(date_time)
        evl_dates_dic[os.path.basename(file)] = (date_time[0], date_time[-1])
    
    date_raw_dic={}
    for file in raw_files_list:
        file = os.path.basename(file)
        cruise, date, time = file.replace(".raw", "").split("-") # Break apart .raw filename to get cruise, date, and time
        start_dt = np.datetime64(datetime.strptime(date[1:] + time[1:], "%Y%m%d%H%M%S")) # Start of acoustic data collection for file
        date_raw_dic[start_dt] = file # Dictionary of start times and .raw files


    evl_raw_dic = {}
    start_delay = np.timedelta64(30, 'm') # Time delay for starting CTD data collection after acoustics 
    for evl_file in evl_dates_dic:
        (stime, etime) = evl_dates_dic[evl_file]
        raw_list = []
        for date in date_raw_dic:
            # If acoustic start is between CDT times or start_delay before, record that raw file
            if (date > stime and date < etime) or abs(date - stime) < start_delay:
                raw_list.append(date_raw_dic[date])
        evl_raw_dic[evl_file] = raw_list
    
    if len(outfile_path) != 0:
        outfile_path = os.path.normpath(outfile_path + "/evl_raw_matches.list")
        outfile = open(outfile_path, 'w')
        headerline = 'Cruise ' + cruise + " CTD and matching acoustic (.raw) files"
        outfile.write(headerline+'\n')
        for evl in evl_raw_dic:
            outfile.write(evl + " " + " ".join(evl_raw_dic[evl]) + "\n") # write CDT file (.evl) followed by overlapping .raw files
        outfile.close()

    return evl_raw_dic

def evl_raw_dic_from_file(evl_raw_file):
    '''
    evl_raw_dic_from_file: reads in file of evl raw matchs and
                           creates a dictionary
    Input: evl_raw_file (string file name and path): path/filename of evl raw match file created by match_raw_evl
           function:
           Format: example_cast.evl raw_overlap_1.raw raw_overlap_2.raw
           Example: ctd017.evl GU19_05-D20191027-T140551.raw GU19_05-D20191027-T143423.raw

           There is one line per CTD cast
    Output: Dictionary with .evl file names as keys with a list of .raw files that overlap in time as the keys
    '''
    evl_raw_file = os.path.normpath(evl_raw_file)
    evl_raw_dic={}
    with open(evl_raw_file, 'r') as file_matches:
        evl_raw = file_matches.readlines()
    for rows in evl_raw[1:]:
        evl, raw = rows.split(maxsplit=1)
        evl_raw_dic[evl] = raw.split()
    return evl_raw_dic
