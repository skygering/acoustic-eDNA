import matplotlib.pyplot as plt
import numpy as np
import json
import math
from math import isclose
from numpy.ma import mod, sort
from scipy.signal import savgol_filter
from echolab2.processing import line
from statistics import mean
import matplotlib.dates as mdates
import os

def read_casts(infile):
    '''
    Creates a dictionary of casts that lists the depths of the eDNA samples
    '''
    with open(infile, 'r') as cast_file:
        cast_data = cast_file.readlines()
    cast_dic = {}
    for cast in cast_data[1:]: # excludes header
        num, depths = cast.split(maxsplit=1)
        depths = [to_float(x) for x in depths.strip().split()]
        cast_dic[int(num)] = depths
    return cast_dic

def calc_slope(x1, y1, x2, y2):
    '''
    Calculated and returns the slope from two x and two y values representing two points
    Assumes points have different x values
    '''
    return (y2 - y1)/(x2 - x1) # what if data is a "bad" data point?

def new_segment(segments, seg_num, bottle):
    ''' 
    Creates a new segment - new entry in segment dictionary
    '''
    segments[seg_num] = {'bottle': bottle, "depth": math.nan, "usable": False, "points" : []}
    return segments

def to_float(x):
    '''
    Attempts to turn value x into a float - returns NaN otherwise
    '''
    try:
        val = float(x)
    except:
        val = math.nan
    return val

def create_segments(atol_zero, x, y):
    x_hat = [num.astype("float") for num in x]
    y_hat = savgol_filter(y, 21, 1)

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


def plot_segments(segments, title = "CTD Track: Depth vs Time"):
    fig, ax = plt.subplots()
    for num in segments:
        style = 'r-'
        if segments[num]["bottle"]:
            style = 'b-'
            if segments[num]["usable"]:
                style = 'b:'
        x_seg, y_seg = zip(*segments[num]["points"])
        x_seg = [np.datetime64(date) for date in x_seg]
        ax.plot(x_seg, y_seg, style)
    fig.autofmt_xdate()
    xfmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_title(title)
    ax.set_xlabel("Time (H:M)")
    ax.set_ylabel("Depth (m)")
    ax.invert_yaxis()
    plt.show()


def usable_samples(segments, atol_depth, casts):
    for num in range(len(segments)):
        seg = segments[num]
        if seg["bottle"]:
            times, depths = zip(*seg["points"])
            seg["depth"] = mean(depths)
            for sample_depth in casts:
                if isclose(seg["depth"], sample_depth, abs_tol=atol_depth) and seg["depth"] >= 5 and sample_depth >=5: # 5s are transducer depth - make not magic number
                    seg["usable"] = True
    return segments

def check_segments(segments, x, y):
    print("Number of segments found with this parameter: ")
    print(len(segments))
    print("The CTD track should be seperated into ascents/descents (red) and plateaus (blue). The plateaus are where data was taken. \n\
Do the number of segments above match the total number of ascents/descents/plateaus and are they tagged with the right color? [y/n]")
    plot_segments(segments)
    correct_graph = input()
    plt.close()
    if correct_graph.lower() != "y":
        print("Please enter a new value for the absolute tolerance. Go up by 0.00001 if too many segments. Go down by same amount if too few.")
        new_atol_zero = float(input())
        segments = create_segments(new_atol_zero, x, y)
        segments = check_segments(segments, x, y)
    return segments

def check_samples(segments, casts):
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
    plot_segments(segments)
    correct_graph = input()
    plt.close()
    if correct_graph.lower() != "y":
        print("Please enter a new value for absolute tolerance to find eDNA sample depths. Go up by 0.5 to 1m if too few samples. Go down by same amount if too many.")
        new_atol_depth = float(input())
        for num in segments:
            segments[num]["usable"] = False
        segments = usable_samples(segments, new_atol_depth, casts)
        segments = check_samples(segments, casts)
    return segments



def make_segments_json(eDNA_cast_depths, ctd_evl_file, outfile_path):
    plt.ion()
    depth_data = line.read_evl(ctd_evl_file)
    x=depth_data.ping_time
    y=depth_data.data
    print("***Analyzing " + os.path.basename(ctd_evl_file).replace(".evl", "") + "***")
    print("SEGMENTATION")
    # adjust for best results
    atol_zero = max(y)*3.16E-7 + 9.6E-5
    print("Absolute Tolerance (maximum absolute distance a slope can be from 0 to be classified as a plateau): ")
    print(round(atol_zero, 5))

    segments = create_segments(atol_zero, x, y)

    # setup for segments

    segments = check_segments(segments, x, y)

    print("\nSAMPLE DEPTHS")
    atol_depth = 2
    print("Absolute tolerance of eDNA sample depths (maximim distance plateu's mean depth differ \
from recorded sample depth to be classified as sample site): ")
    print(atol_depth)
    
    segments = usable_samples(segments, atol_depth, eDNA_cast_depths)
    segments = check_samples(segments, eDNA_cast_depths)
    
    print("Saving .json file with segment classification.\n")
    with open(outfile_path, 'w') as outfile: # might be better to only save the  "usable" segments and then note the min/max depth and min/max time
        json.dump(segments, outfile)


