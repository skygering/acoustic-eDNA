import matplotlib.pyplot as plt
import numpy as np
import json
import math
from math import isclose
from numpy.ma import mod
from scipy.signal import savgol_filter
from echolab2.processing import line
from statistics import mean

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

casts = cast_dic[3] # which cast
ctd = '/Volumes/GeringSSD/GU201905_CTD/ctd003.evl' # which CTD
outfile_path = '/Volumes/GeringSSD/segments_cdt003.json'
depth_data = line.read_evl(ctd)

x=depth_data.ping_time
y=depth_data.data
x_hat = [num.astype("float") for num in x]
window_len = int(len(x) * 0.05)

y_hat = savgol_filter(y, 21, 1)


# adjust for best results
atol_zero = max(y)*3.16E-7 + 9.6E-5
print("Absolute Total (atol):")
print(round(atol_zero, 5))

# setup for segments
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

print("Number of segments with these parameters:")
print(len(segments))

print("Close out of the graph once you have determined if it has the correct number of segments.")
plot_segments(segments)

print("Did the graph look correct?")
input1 = input()
print(input1)

atol_depth = 2
print("Depth variance to find eDNA casts:")
print(atol_depth)

for num in range(len(segments)):
    seg = segments[num]
    if seg["bottle"]:
        times, depths = zip(*seg["points"])
        seg["depth"] = mean(depths)
        for cast_depth in casts:
            if isclose(seg["depth"], cast_depth, abs_tol=atol_depth) and cast_depth >= 5:
                seg["usable"] = True



plot_segments(segments)


with open(outfile_path, 'w') as outfile:
    json.dump(segments, outfile)