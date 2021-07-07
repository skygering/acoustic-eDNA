from echolab2.processing import line
from datetime import date, datetime, timedelta
from statistics import mean
import numpy as np
import os
import math
import glob
from math import isclose
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from numpy.ma import mod, sort
import json
import CTD_EK_plotting as plotting

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

def to_float(x):
        '''
        to_float: Attempts to turn value x into a float - returns NaN otherwise
        '''
        try:
            val = float(x)
        except:
            val = math.nan
        return val

def read_casts(infile):
    '''
    read_casts: reads in a file that notes the cast and the depths at which usable samples (eDNA samples) were taken
    Input: infile (string) - path + filename to a file with cast and depths of eDNA samples
                             Format: cast_num eDNA_depth
                             Example: 1 10 47
                             Note: There is a header ("cast eDNA")
    Output: cast_dic (dictionary) - dictionary with keys as cast numbers and values as a list of depths were eDNA samples were taken
    '''
    with open(infile, 'r') as cast_file:
        cast_data = cast_file.readlines()
    cast_dic = {}
    for cast in cast_data[1:]: # excludes header
        num, depths = cast.split(maxsplit=1)
        depths = [to_float(x) for x in depths.strip().split()]
        cast_dic[int(num)] = depths
    return cast_dic


def create_segments_dic(x, y, atol_zero):
    '''
    create_segments_dic: takes depth (y) over time (x) data from a CTD trace and splits it into segments by slope to
                         identify plateau's in the trace where water bottle data was collected
    Inputs: x (list of numpy64 datetime objects) - times in the CTD trace
            y (list of floats) - depths of the CTD trace over time
            atol_zero (float) - magnitude of the the minimum absolute tolerance when determining if two numbers are approimatly zero -
                                number will probably be suprisingly small - try 0.00002 to start
    Outputs: segment dictionary - nested dictionary; outermost keys are segment number (i.e. 0, 1, 2), followed by the following structure for each
                                  segment:
                                  Format: 'bottle': True/False, "depth": float/math.nan, "usable": False, "points" : [all points]

                                  -If the slope of a segment is approximatly 0, as determined by the atol_zero parameter, then a water bottle 
                                  sample was collected so 'bottle' is True, else False
                                  -If 'bottle' is True, the mean depth of the segment is recorded under 'depth' key, else it is math.nan
                                  -The 'usable' keys are all set to false here - they can be set to mark eDNA sample locations with the 
                                   mark_usable_segments function below
                                  -The 'points' key is followed by tuples of time/depth points (x,y) that are in the segment
    '''

    def calc_slope(x1, y1, x2, y2):
        '''
        calc_slope: calculate and return the slope from two x and two y values representing two points
        Inputs: x1 (float/int): x value of first point; y1 (float/int) - y value of first point
                x2 (float/int): x value of second point; y2 (float/int) - y value of second point
        Outputs: float slope between the two points
        Warning: Assumers that the x values of the two points are different to avoid divide by 0 error
        '''
        return (y2 - y1)/(x2 - x1)

    def new_segment(segments, seg_num, bottle):
        ''' 
        new_segment: adds a new segment to the segment dictionary
        Inputs: segments (nested dictionary of segments) - current dictionary
                seg_num (int) - number of new segment; will be the key for new segment
                bottle (boolean) - True if slope of new segment is approximatly 0, False otherwise
        Outputs: segments - updated segment dictionary with new added segment
        '''
        segments[seg_num] = {'bottle': bottle, "depth": math.nan, "usable": False, "points" : []}
        return segments
    # create_segments_dic starts here
    x_hat = [num.astype("float") for num in x] # turns numpy64 datetime objects into floats
    y_hat = savgol_filter(y, 21, 1) # smooths the data into somewhat linear segments

    seg_num = 0
    segments = {}
    current_slope = calc_slope(x_hat[0], y_hat[0], x_hat[1], y_hat[1])
    segments = new_segment(segments, seg_num, isclose(current_slope, 0, abs_tol= atol_zero))

    for j in range(len(x)-1):
        x_j = np.datetime_as_string(x[j])
        y_j = float(y[j])
        new_slope = calc_slope(x_hat[j], y_hat[j], x_hat[j+1], y_hat[j+1])
        new_bottle = isclose(new_slope, 0, abs_tol= atol_zero)
        if segments[seg_num]["bottle"] == new_bottle: # we are still on the same segment - be that plateau or otherwise
            segments[seg_num]["points"].append((x_j, y_j)) # add in current point
        else: # we have now switched to a new type of segment
            if len(segments[seg_num]["points"]) > 5: # if the previous segment is more than 2 points long
                seg_num += 1
                segments = new_segment(segments, seg_num, new_bottle) # makes empty new segment
                segments[seg_num]["points"].append((x_j, y_j)) # Adds current point to new segment
            else: # short previous segment - due to noise (too short to be actual segment)
                segments[seg_num-1]["points"] += segments[seg_num]["points"] # add glitch points to previous segment
                segments[seg_num-1]["points"].append((x_j, y_j)) # add new point to previous segment 
                # since we only switch segments when they are diff type, previous and current must be same type with "glitch" inbetween 
                del segments[seg_num] # delete "glitch segment"
                seg_num -= 1
    segments[seg_num]["points"].append((np.datetime_as_string(x[-1]), float(y[-1]))) # add the last point to last segment
    return segments

def mark_usable_depth(segments, cast_depths, transducer_depth, atol_depth = 2):
    '''
    mark_usable_segments: given a segment dictionary and a list of depths of eDNA samples for each cast, mark the 'usable' tag in the
                          segment dictionary True if a segment is at an eDNA sample depth and leave it as false otherwise
    Inputs: segments (nested dictionary) - segment dictionary created by create_segments_dic function
            casts (numeric list) - list of depths eDNA samples were taken at during current cast
            atol_depth (float) - the minimum absolute difference in meters that the mean depth of the segment can differ from the depth
                                 specified in the cast file
    Outputs: segment dictionary with "usable" tags updated and "depth" 
    '''
    for num in range(len(segments)):
        seg = segments[num]
        if seg["bottle"]: # only segments with water samples taken could be a eDNA sample site
            times, depths = zip(*seg["points"])
            seg["depth"] = mean(depths) # calculates mean depth of segments and saves in segments
            for sample_depth in cast_depths:
                if isclose(seg["depth"], sample_depth, abs_tol=atol_depth) and seg["depth"] >= transducer_depth and sample_depth >=transducer_depth:
                    seg["usable"] = True # sets "usable" tag based on sample depth
    return segments


def interactive_segment_maker(eDNA_cast_depths, ctd_evl_file, transducer_depth, outfile_path=""):
    '''
    interactive_segment_maker: interactive script that walks user through creating .json of segments (also denotes which are usable)
                               for a single cast - if minimum absolute tolerance variables are wrong, script allows you to adjust and 
                               dynamically changes dictionary before saving
    Inputs: eDNA_cast_depths (integer list) - list of depths eDNA data was collected at for cast
            ctd_evl_file (string) - filename and path to CTD trace .evl file
            transducer_depth (float) - depth of transducer in meters
            outfile_path (string) - optional filename and path to save resultant .json file to
    Outputs: segment dictionary - dictionary for a single cast - outermost key is segment number (i.e. 0, 1, 2) while inner keys are:
                                    Format: 'bottle': True/False, "depth": float/math.nan, "usable": False, "points" : [all points]
            If outfile_path is provided, the .json file is saved to the outfile path with the name of the ctd_evl_file with "_segments"
            appended to the end of the name and the .json file extension
    Note: This is RECOMMENDED for making .json files - it is much easier than doing all of them with create_segments_dic
    Warning: This depends on a formula created to estimate good values for minimum absolute tolerance (atol_zero):
                atol_zero = max(y)*3.16E-7 + 9.6E-5
             This might need to be updated for new datasets - manually find good values for atol_zero using create_segments_dic by trying 
             values for a few casts -  make a graph of atol_zero vs maximum depth of cast and use the fit line
    '''

    def check_segments_dic(segments, x, y):
        '''
        check_segments_dic: interactive check of create_segments_dic function - allows user to confirm that the correct number of segments 
                            were found and try a new minimum absolute tolerance if not
        Input: segments (segment dictionary) - segment dictionary after running create_segments_dic()
               x (numpy datetime64 list) - datetime objects from original CTD trace
               y (float list) - depth values from original CTD trace
        Output: allows dynamic adjustment of minimum absolute tolerance - once user is satisfied updated segment dictionary is returned
        '''
        # check on number of segments
        print("Number of segments found with this parameter: ")
        print(len(segments))
        print("The CTD track should be seperated into ascents/descents (red) and plateaus (blue). The plateaus are where data was taken. \n\
Do the number of segments above match the total number of ascents/descents/plateaus and are they tagged with the right color? [y/n]")
        plotting.plot_segments(segments)
        correct_graph = input() # check if all segments are found
        plt.close()
        if correct_graph.lower() != "y": # all segments were not found
            print("Please enter a new value for the absolute tolerance. Go up by 0.00001 if too many segments. Go down by same amount if too few.")
            new_atol_zero = float(input())
            segments = create_segments_dic(x, y, new_atol_zero) # try to find segments again
            segments = check_segments_dic(segments, x, y) # check new segments
        return segments

    def check_usable_segments(segments, casts):
        '''
        check_usable_segments: interactive check of mark_usable_depth function - allows user to confirm that the depth of eDNA samples
                               were found and try a new minimum absolute tolerance for depth if not
        Inputs: segments (segment dictionary) - segment dictionary after running mark_usable_depth()
                casts (integer list) - list of depths at which eDNA samples were taken
        Output: allows dynamic adjustment of minimum absolute tolerance to find eDNA samples segments
                - once user is satisfied updated segment dictionary is returned
        '''
        # check if all eDNA samples are identified
        print("Expected, usable (below transducer depth) eDNA sample depths: ")
        print([c for c in casts if c > 5]) # this 5 is transducer depth - make not a magic number
        found_samples = []
        for num in segments:
            seg = segments[num]
            if seg["usable"]:
                found_samples.append(round(seg["depth"], 2))
        print("Found, usable eDNA sample depths in segments: ")
        print((sort(found_samples)))

        print("The usable eDNA samples found should be displayed above and identifed with dotted lines on the graph. Are all expected eDNA sample depths found? [y/n]")
        plotting.plot_segments(segments)
        correct_graph = input() # check if all samples are identified
        plt.close()
        if correct_graph.lower() != "y": # if all samples are not identified
            print("Please enter a new value for absolute tolerance to find eDNA sample depths. Go up by 0.5 to 1m if too few samples. Go down by same amount if too many.")
            new_atol_depth = float(input())
            for num in segments:
                segments[num]["usable"] = False # reset all "usable" segments to unusable before trying to find them again
            segments = mark_usable_depth(segments, casts, transducer_depth, new_atol_depth) # try to find samples
            segments = check_usable_segments(segments, casts) # check samples
        return segments

    # interactive_segment_maker starts here
    plt.ion()
    depth_data = line.read_evl(ctd_evl_file)
    x=depth_data.ping_time
    y=depth_data.data

    print("***Analyzing " + os.path.basename(ctd_evl_file).replace(".evl", "") + "***")
    print("SEGMENTATION")
    atol_zero = max(y)*3.16E-7 + 9.6E-5 # might need to adjust for best results
    print("Absolute Tolerance (maximum absolute distance a slope can be from 0 to be classified as a plateau): ")
    print(round(atol_zero, 5))

    segments = create_segments_dic(x, y, atol_zero) # find segments
    segments = check_segments_dic(segments, x, y) # check segments

    print("\nSAMPLE DEPTHS")
    atol_depth = 2
    print("Absolute tolerance of eDNA sample depths (maximim distance plateau's mean depth differ \
from recorded sample depth to be classified as sample site): ")
    print(atol_depth)
    
    segments = mark_usable_depth(segments, eDNA_cast_depths, transducer_depth, atol_depth) # find eDNA samples
    segments = check_usable_segments(segments, eDNA_cast_depths) # check samples
    
    if len(outfile_path) != 0: # if we have an outfile, save outfile
        print("Saving .json file with segment classification.\n")
        with open(outfile_path, 'w') as outfile: # might be better to only save the  "usable" segments and then note the min/max depth and min/max time
            json.dump(segments, outfile)

    return segments





