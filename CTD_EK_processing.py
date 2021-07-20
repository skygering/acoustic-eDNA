from echolab2.plotting.matplotlib import echogram
from echolab2.processing import line, processed_data
from echolab2.instruments import EK80
from datetime import datetime, timedelta
from statistics import mean
import numpy as np
import os
import math
from math import isclose
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from numpy.ma import sort
import json
import CTD_EK_plotting as plotting
import glob
from itertools import combinations
import warnings


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
                             Note: There is NO header
    Output: cast_dic (dictionary) - dictionary with keys as cast numbers and values as a list of depths were eDNA samples were taken
    '''
    with open(infile, 'r') as cast_file:
        cast_data = cast_file.readlines()
    cast_dic = {}
    for cast in cast_data:
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
            If outfile_path is provided, the .json file is saved to the outfile path
    Note: This is RECOMMENDED for making .json files - it is much easier than doing all of them with create_segments_dic
    Warning: This depends on a formula created to estimate good values for minimum absolute tolerance (atol_zero):
                atol_zero = max(y)*3.16E-7 + 9.6E-5
             This might need to be updated for new datasets - manually find good values for atol_zero using create_segments_dic by trying 
             values for a few casts -  make a graph of atol_zero vs maximum depth of cast and use the fit line
    '''
    def get_atol(atol):
        '''
        allows user to try multiple times to give new atol values and is used for error handling with bad entries
        '''
        try: 
            atol = atol.split()
            if len(atol) != 1:
                raise
            new_atol =  float(atol[0])
        except:
            new_atol = input("Did not enter a valid offset. Please try again:")
            new_atol =  get_atol(new_atol)
        return new_atol

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
        plotting.plot_segments(segments)
        correct_graph = input("The CTD track should be seperated into ascents/descents (red) and plateaus (blue). The plateaus are where data was taken. \n\
Do the number of segments above match the total number of ascents/descents/plateaus and are they tagged with the right color? [y/n]") # check if all segments are found
        plt.close()
        if correct_graph.lower() != "y": # all segments were not found
            new_atol_zero = get_atol(input("Please enter a new value for the absolute tolerance. Go up by 0.00001 if too many segments. Go down by same amount if too few."))
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

        plotting.plot_segments(segments)
        correct_graph = input("The usable eDNA samples found should be displayed above and identifed with dotted lines on the graph. Are all expected eDNA sample depths found? [y/n]") # check if all samples are identified
        plt.close()
        if correct_graph.lower() != "y": # if all samples are not identified
            new_atol_depth = get_atol(input("Please enter a new value for absolute tolerance to find eDNA sample depths. Go up by 0.5 to 1m if too few samples. Go down by same amount if too many."))
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
        with open(outfile_path, 'w') as outfile:
            json.dump(segments, outfile)

    return segments


def raw_to_Sv(ek, fq, transducer_offset):
    '''
    raw_to_Sv: function takes an ek object and returns the Sv processed data object and
               the calibration object for the specified frequency and transducer offset
    Inputs: ek (EK80 data object) - see Rick Towler's pyEcholab code for more details
            fq (integer) - frequency for desired data
            transducer_offset (double) - offset of transducer from water surface in meters
    Outputs: Sv object - Sv data object for given frequency and EK object
             Calibration object - calibration object with given transducer offset set
    Note: for more details on Sv and Calibration objects see Rick Towler's pyEcholab code
    '''
    raw_data_list = ek.get_channel_data(frequencies=fq)
    raw_data = raw_data_list[fq][0] 
    for file in raw_data_list[fq][1:]:
        raw_data.append(file)
    cal_obj = raw_data.get_calibration()
    cal_obj.transducer_offset_z = transducer_offset
    return raw_data.get_Sv(calibration=cal_obj, return_depth=True), cal_obj


def crop_Sv(Sv, xlim, ylim):
        '''
        crop_Sv: given x and y index limits, subset an Sv processed data object from pyEcholab
        Inputs: Sv (processed data object from pyEcholab) - run get_Sv (pyEcholab) or raw_to_Sv (see above)
                xlim (integer list) - [lower limit x index, higher limit x index] 
                                      subsets inclusive of both of these limits 
                ylim (integer list) - [lower limit y index, higher limit y index] 
                                      subsets inclusive of both of these limits 
        Outputs: cropped Sv data object with all major data attributes (data, ping_time, deoth, 
                 transducer_offset, n_pings, heave) cropped to given limits
        '''
        # data has dims (ping_time, depth) - rows are ping_time, columns are depth
        Sv.data = Sv.data[xlim[0]:xlim[1]+1, ylim[0]:ylim[1]+1] 
        Sv.ping_time = Sv.ping_time[xlim[0]:xlim[1]+1]
        Sv.depth = Sv.depth[ylim[0]:ylim[1]+1]
        Sv.transducer_offset = Sv.transducer_offset[xlim[0]:xlim[1]+1]
        Sv.n_pings = len(Sv.ping_time)
        Sv.heave = Sv.heave[xlim[0]:xlim[1]+1]
        return Sv # is there a better way to resize this so that everything gets resized?

def subset_segments_Sv(segment_dic, raw_files, toffsets, doffsets, transducer_offset, frequencies = [18000, 38000, 120000, 200000]):
    '''
    subset_segments_Sv: takes a segment dictionary generated by create_segments_dic - for each usable segment, the matching
                        raw file is subset around the segment using time and depth offsets to create a smaller rectangle of
                        the raw data for each frequency
    Inputs: segment_dic (segment dictionary) - segment dictionary after running create_segments_dic or interactive_segment_maker
            raw_files (list of string) - list of filenames of raw files that correspond to CTD profile in segment_dic
            toffsets (integer list) - two item list with minute offsets from the start of data collection in segments 
                                      first item for minutes before data collection and second item for minutes after start of collection
            doffsets (double list) - two item list with meter offsets from mean depth of segment
                                     first item is meters above mean, second is meters below mean
            transducer_offset (double) - offset of transducer from water surface in meters
            frequencies (integer list) - list of frequnecies within raw data
    Outputs: (1) nested dictionary for each usable segment with points within defined subset for each frequency
             Format: 
                        {depth 1: {frequency 1: cropped Sv data object for frequency 1 at depth 1
                                   frequency 2: cropped Sv data object for frequency 2 at depth 1
                                   frequency 3: cropped Sv data object for frequency 3 at depth 1
                                   frequency 4: cropped Sv data object for frequency 4 at depth 1
                                  }
                         depth 2: {.
                                   .
                                   .}
                         ...
                        }
    Note: This is for one CTD profile - if you want to do more than one you need to loop
            '''
    ek80 = EK80.EK80()
    ek80.read_raw(raw_files) # raw files that will be subset
    subset_segment_dic = {}
    for num in segment_dic:
        seg = segment_dic[num]
        if seg["usable"]: # we only want to subset depths at which eDNA data was taken
            y_mean = float(seg["depth"]) # depth of sample
            fq_dic = {}
            subset_segment_dic[round(y_mean)] = fq_dic # outermost keys are sample depth
            points_time, points_depth = zip(*seg["points"])
            points_time = [np.datetime64(x) for x in points_time]

            # bounds of subset
            y_min = y_mean - doffsets[0]
            y_max = y_mean + doffsets[1]
            x_min = min(points_time) - np.timedelta64(toffsets[0], 'm') 
            x_max = min(points_time) + np.timedelta64(toffsets[1], 'm')

            Sv18, cal18 = raw_to_Sv(ek80, 18000, transducer_offset)
            # find indices closest to bounds
            # once we resample, these should be the same for each frequency
            depths = np.array(Sv18.depth)
            y_min_idx = np.argmin(np.abs(depths - y_min))
            y_max_idx = np.argmin(np.abs(depths - y_max))

            times = np.array(Sv18.ping_time)
            x_min_idx = np.argmin(np.abs(times - x_min))
            x_max_idx = np.argmin(np.abs(times - x_max))

            for fq in frequencies:
                Sv, cal =  raw_to_Sv(ek80, fq, transducer_offset) # get Sv data for each frequency
                if cal.sample_interval != cal18.sample_interval: # resample to match 18kHz
                    Sv.match_samples(Sv18)
                subset_Sv = crop_Sv(Sv.copy(), (x_min_idx, x_max_idx), (y_min_idx, y_max_idx)) # "crop" Sv to get subset
                fq_dic[fq] = subset_Sv # Sv object saved

    return subset_segment_dic

def subset_to_json(subset_dic):
    '''
    subset_to_json: takes  a subset dictionary created by subset_segments_Sv and turns it into a format that can be dumped to a .json file
    Inputs: subset_dic (subset dictionary of Sv objects)
    Outputs: (1) Dictionary with same structure as subset_dic but instead of the innermost values being Sv data objects from pyEcholab,
             they're nested lists of doubles - this is just the Sv data and does not have the ping_time or depth data like the Sv object
             (2) Nested dictionary with depths for outermost keys and "depths" and "ping time" for inner keys with the depth and ping time
                 values for each list stored as 1D arrays 
    '''
    bounds_dic = {}
    for depth in subset_dic:
        depth_subset = subset_dic[depth]
        frequencies = list(depth_subset.keys())
        pings_arr = [np.datetime_as_string(x) for x in depth_subset[frequencies[0]].ping_time]
        depth_arr = depth_subset[frequencies[0]].depth.tolist()
        bounds_dic[depth] = {"depth": depth_arr, "ping_time": pings_arr}
        for fq in frequencies:
            depth_subset[fq] = depth_subset[fq].data.tolist() # lists can be "dumped" to .json files
    return subset_dic, bounds_dic


def interactive_subset_maker(segment_dic, raw_files, transducer_offset, toffsets, doffsets,
                             frequencies = [18000, 38000, 120000, 200000], thresholds = [-90, -20], outfile_path=""):
    '''
    interactive_subset_maker: interactive check of subset_segments_Sv function - allows user to confirm the bounds of the subsets for each depth
                              and to choose to remove any frequencies with bad data
    Inputs: segment_dic (segment dictionary) - segment dictionary after running create_segments_dic or interactive_segment_maker
            raw_files (list of string) - list of filenames of raw files that correspond to CTD profile in segment_dic
            toffsets (integer list) - two item list with minute offsets from the start of data collection in segments 
                                      first item for minutes before data collection and second item for minutes after start of collection
            doffsets (double list) - two item list with meter offsets from mean depth of segment
                                     first item is meters above mean, second is meters below mean
            transducer_offset (double) - offset of transducer from water surface in meters
            frequencies (integer list) - list of frequnecies within raw data
            thresholds (integer list) - thresholds for echogram plotting with pyEcholab
            outfile_path (string) - if a path and filename is provided, a json file is saved that has Sv data points in an array (not Sv object)
    Outputs: nested dictionary for each usable segment with points within defined subset for each frequency - see subset_segments_Sv function
             comments for format
    Note: This is for one CTD profile - if you want to do more than one you need to loop
    '''


    def check_subset_bounds(subset, depth_val):
        '''
        check_subset_bounds: checks the bounds of subset for one depth 
        Inputs: subset (dictionary) - dictionary of ONE depth, which has keys that are frequencies and values that are subset Sv objects
                depth_val - depth of the subset
        Output: dictionary of specific depth with bounds specified and approved by user
        '''
        def get_offset(message, type):
            '''
            get_offset: allows user to try multiple times to give new offsets and is used for error handling with bad entries
            '''
            offsets = []
            try:
                offsets = input(message).split(",")
                if len(offsets) !=2:
                    raise
                if type == "int":
                    offsets = [int(x) for x in offsets]
                else: offsets = [float(x) for x in offsets]
            except: 
                print("You input was wrong. Wrong number or type of inputs. Time offsets must be integers, depth can be decimals. Please try again")
                offsets = get_offset(message, type)
            return offsets


        fig, ax = plt.subplots(constrained_layout = True)
        ax.set_title("Example Subset at " + str(depth_val) + "m")
        echogram.Echogram(ax, subset[frequencies[0]], threshold = thresholds)
        correct_bounds = input("Are the bounds of the subset postioned well (i.e does not go into the surface, the bottom)? [y/n]")
        plt.close()
        if correct_bounds.lower() != "y":
            new_toffsets = get_offset("Time offsets in minutes (minutes before and after start of data capture; seperate with a comma):", "int")
            new_doffsets = get_offset("Depth offsets in meters (meters above and below CTD sample depth; seperate with a comma):", "float")
            segment = {}
            n_segment = 0
            for seg_num in segment_dic:
                seg = segment_dic[seg_num]
                if seg["usable"] and round(seg["depth"]) == int(depth_val):
                    segment = seg
                    n_segment = seg_num
            subset= subset_segments_Sv({n_segment: segment}, raw_files, new_toffsets, new_doffsets, transducer_offset)
            subset = check_subset_bounds(subset[depth_val], depth_val)
        return subset
    
    def check_subset_freq(subset, depth_val):
        '''
        check_subset_freq: checks if all frequencies should be included and gives users the chance to exclude 
                           frequencies from subset dictionary
        Inputs: subset (dictionary) - dictionary of ONE depth, which has keys that are frequencies and values that are subset Sv objects
                depth_val - depth of the subset
        Output: dictionary of specific depth with all frequencies remaining; freqencies specified by user removed
        '''
        def remove_fq(subset):
            '''
            remove_fq: allows users to remove frequencies from dictionary
            '''
            rm_fq = input("List freqencies you want to remove seperated by commas: ").split(",")
            for fq in rm_fq:
                try:
                    subset.pop(int(fq)) # remove specified key
                except: # user did not specify a key within the dictionary
                    again = input(str(fq) + " is not a frequency. It could not be removed. Do you want to try another frequency? [y/n]")
                    if again.lower() == "y":
                        remove_fq(subset) # user can try again
            return subset

        n_fq = len(frequencies)
        fig, ax = plt.subplots(n_fq, figsize=(10,8), constrained_layout = True)
        fig.suptitle("Subset at " + str(depth_val))
        for f in range(n_fq): # loop through freqencies
            echogram.Echogram(ax[f], subset[frequencies[f]], threshold = thresholds) # plot to see all freqencies
            ax[f].set_title(str(frequencies[f]) + "Hz")
        all_fq = input("Do you want all frequencies shown included in the output dictionary? [y/n]")
        if all_fq.lower() != "y":
            subset = remove_fq(subset) # remove frequencies specified by user
        plt.close()
        return subset

    # interactive_subset_maker starts here
    plt.ion()
    subset_dic = subset_segments_Sv(segment_dic, raw_files, toffsets, doffsets, transducer_offset, frequencies) # starting subsets
    eDNA_depths = list(subset_dic.keys())
    n_depths = len(eDNA_depths)
    for d in range(n_depths):
        depth = eDNA_depths[d]
        print("Subsetting at " + str(depth) + "m")
        subset = check_subset_bounds(subset_dic[depth], depth) # check each depth's subset bounds
        subset = check_subset_freq(subset, depth) # check if all freqencies have good data
        subset_dic[depth] = subset # if subset was updated, replace within the subset dictionary

    if len(outfile_path) != 0: # if we have an outfile, save outfile
        print("Saving .json file with subsets (nested lists of data saved rather than Sv objects).")
        subset_dic = subset_to_json(subset_dic) # change subset dictionary into format that can be dumped to .json file
        with open(outfile_path, 'w') as outfile:
            json.dump(subset_dic, outfile)
    return subset_dic


def calc_MFI(Sv_data_dic, delta = 40, bad_fq = [], global_norm = False):
    '''
    calc_MFI: creates processed data object with MFI classification from a dictionary of Sv data with frequencies as keys
    Inputs: Sv_data_dic (Sv data object dictionary) - dictionary with freqencies as keys - can be created with subset_segments_Sv()
            delta (integer) - parameter needed for MFI calculation - see Trenkel, Verena M., and Laurent Berger. "A fisheries acoustic 
                              multi-frequency indicator to inform on large scale spatial patterns of aquatic pelagic ecosystems."
            local_norm (boolean) - calculating MFI requires normalizing the Sv data by maximum value - if True each frequency is maximized
                                   by individual maximum rather than global (between all frequencies) maximum
            bad_fq (integer list) - list of freqencies to NOT use in MFI calculation (i.e. 200000)
    Output: MFI processed data object where the data attribute is the MFI calculation, while the ping_time, n_pings, and depth are the same
            as in input Sv objects
    '''
    f=list(Sv_data_dic.keys()) # freqencies of Sv data availible
    for b in bad_fq:
        try: f.remove(b)
        except: print(str(b) + " was not in frequency list and can't be removed")
    
    f = [int(i)/1000 for i in f]
    f.sort()
    f_inv = {i: 1/i for i in f} # inverse freqency calues
    nf = len(f)

    if nf < 3:
        warnings.warn("Two or less frequencies provided - cannot calcualte MFI!")

    f_comb = list(combinations(f, 2)) # combinations of frequencies
    dist_dic = {}
    for comb in f_comb:
        dist_dic[comb] = 1-math.exp((-1*abs(comb[0]-comb[1]))/delta) # distance calculation
    print(Sv_data_dic.keys())
    sv_data_dic={i: None for i in f} # linearize Sv data
    for i in f:
        if Sv_data_dic[i*1000] is not None:
            sv_data_dic[i] = Sv_data_dic[i*1000].copy()
            sv_data_dic[i].data = 10**(sv_data_dic[i].data/10)
   
    max_vals = np.array([np.amax(sv.data) for sv in sv_data_dic.values()]) # maximum sv value per each freqencys' sv data
    min_vals = np.array([np.amin(sv.data) for sv in sv_data_dic.values()]) # minimum sv value per each freqencys' sv data
    if global_norm: 
        max_vals = np.full(max_vals.shape, np.max(max_vals))
        min_vals = np.full(min_vals.shape, np.min(min_vals))
    print(max_vals)
    
    Norms={i: None for i in f} # normalized sv data for each frequency
    for i in range(len(f)):
        norm_sv = sv_data_dic[f[i]].copy()
        # (sV(fi)-min(sV(f)))/(max(sV(f)) -min(sV(f) normalization equation as confirmed by Berger and Trenkel
        norm_sv.data = ((sv_data_dic[f[i]].data - min_vals[i])/(max_vals[i] - min_vals[i])).astype(float)
        Norms[f[i]] = norm_sv

    # MFI calculatin
    a=[]
    b=[]
    for fq_pair in dist_dic:
        f0 = fq_pair[0]
        f1 = fq_pair[1]
        if Norms[f0] is not None and Norms[f1] is not None:
            x = Norms[f0].data * Norms[f1].data * f_inv[f0] * f_inv[f1]
            a.append(dist_dic[fq_pair] * x)
            b.append(x)
    MFI_data = ((sum(a)/sum(b))-0.4)/0.6 

    # create processed data object from pyEcholab
    MFI_obj = processed_data.processed_data("None", math.nan, "MFI")
    MFI_obj.data = MFI_data
    sv_data = sv_data_dic[f[0]]
    MFI_obj.n_pings = sv_data.n_pings
    MFI_obj.ping_time = sv_data.ping_time
    MFI_obj.depth = sv_data.depth
    return MFI_obj

def processed_data_from_dic(data, bounds_dic, type = "Sv"):
    '''
    processed_data_from_dic: takes a data array and a bound dictionary and creates a processed data object from
                             Rick Towler's pyEcholab code
    Inputs: data (2D data array or numpy array) - usually Sv data
            bounds_dic (dictionary) - bounds for one cast i.e. {"ping_time": array, "depth": array}
            type (string): name for type of data
    Outputs: processed data object with attributes data, ping_time, n_pings, and depth
    '''
    obj = processed_data.processed_data("None", math.nan, type)
    obj.data = np.array(data)
    pings = bounds_dic["ping_time"]
    obj.ping_time = np.array([np.datetime64(x) for x in pings])
    obj.n_pings = len(pings)
    obj.depth = np.array(bounds_dic["depth"])
    return obj

def mask_mfi(mfi, Sv, range):
    '''
    mask_mfi: takes a MFI data array, and a Sv data array and applys a mask where the MFI data is in the range
              and applies the mask to the Sv data
    Inputs: mfi (2D float array) - MFI data array
            Sv (2D floar array) - Sv data array of the same size as mfi
            range (float list) - two item list specifying the upper and lower limits of the range of values that
                                 denote desirable areas of the Sv data
    Outputs: outputs a 2D float array of the Sv data with a mfi mask applied
    '''
    mfi = np.array(mfi)
    mask = np.where((range[0] < mfi) & (mfi < range[1]), 1, 0)
    mask_Sv = np.multiply(mask, np.array(Sv))
    mask_Sv[mask_Sv == 0] = -999
    return mask_Sv

def calc_ABC(Sv_obj):
    '''
    calc_ABC: Calculates the Area Backscattering Coefficent (ABC) for a Sv data object
    Inputs: Sv_obj (pyEcholab processed data object) - Sv processed data object with data and depth attributes
    Outputs: returns ABC value (float) for the Sv data provided 
    '''
    Sv_data = np.array(Sv_obj.data)
    n_points = Sv_data.size
    sv_mean = np.sum(10**(Sv_obj.data/10)) # using all data points
    depths = Sv_obj.depth
    bin_thickness = (max(depths) - min(depths))/len(depths) # does not allow for variable ping depths
    return sv_mean * bin_thickness # return ABC value

