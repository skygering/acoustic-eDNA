from itertools import repeat
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import json
import math
from math import isclose
from scipy.linalg.lapack import _check_work_float
from scipy.signal import savgol_filter
from echolab2.processing import line
from statistics import mean

def read_casts(infile):
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
    return (y2 - y1)/(x2 - x1)

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

def plot_segments(segments):
    for num in segments:
        style = 'r-'
        if segments[num]["bottle"]:
            style = 'b-'
            if segments[num]["usable"]:
                style = 'b:'

        x_seg, y_seg = zip(*segments[num]["points"])
        x_seg = [np.datetime64(date) for date in x_seg]
        plt.plot(x_seg, y_seg, style)
    plt.show()
    plt.close()

# Read in depths at which eDNA measurments were taken for each cast
cast_dic = read_casts('/Volumes/GeringSSD/eDNA_cast.txt')


casts = cast_dic[17] # which cast
ctd = '/Volumes/GeringSSD/GU201905_CTD/ctd017.evl' # which CTD
outfile_path = '/Volumes/GeringSSD/segments_cdt017.json'
depth_data = line.read_evl(ctd)
x=depth_data.ping_time
y=depth_data.data


x_hat = [num.astype("float") for num in x]
y_hat = savgol_filter(y, 17, 1)

# adjust for best results
atol = .0003

# setup for segments
seg_num = 0
segments = {}
current_slope = calc_slope(x_hat[0], y_hat[0], x_hat[1], y_hat[1])
segments = new_segment(segments, seg_num, isclose(current_slope, 0, abs_tol= atol))

for j in range(len(x)-1):
    x_j = np.datetime_as_string(x[j])
    y_j = float(y[j])
    new_slope = calc_slope(x_hat[j], y_hat[j], x_hat[j+1], y_hat[j+1])
    new_bottle = isclose(new_slope, 0, abs_tol= atol)
    if segments[seg_num]["bottle"] == new_bottle:
        segments[seg_num]["points"].append((x_j, y_j))
    else:
        seg_num += 1
        segments = new_segment(segments, seg_num, new_bottle) # makes empty new segment
        segments[seg_num]["points"].append((x_j, y_j)) # new segment with point
        current_slope = new_slope # approximate slope of new segment
segments[seg_num]["points"].append((np.datetime_as_string(x[-1]), float(y[-1]))) # add the last point to last segment


for num in range(len(segments)):
    seg = segments[num]
    if seg["bottle"]:
        times, depths = zip(*seg["points"])
        seg["depth"] = mean(depths)
        for cast_depth in casts:
            if isclose(seg["depth"], cast_depth, abs_tol=2) and cast_depth >= 5:
                seg["usable"] = True 

# Plotting
print("Number of segments:")
print(len(segments))
plot_segments(segments)
print("If number of segments is too high, try lowering atol.")

with open(outfile_path, 'w') as outfile:
    json.dump(segments, outfile)