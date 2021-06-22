from itertools import repeat
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import json
import math
from math import isclose
from scipy.linalg.lapack import _check_work_float
from scipy.signal import savgol_filter
from echolab2.processing import line, processed_data
from statistics import mean

def calc_slope(x1, y1, x2, y2):
    '''
    Calculated and returns the slope from two x and two y values representing two points
    Assumes points have different x values
    '''
    return (y2 - y1)/(x2 - x1)

def new_segment(segments, seg_num, slope):
    ''' 
    Creates a new segment - new entry in segment dictionary
    '''
    ctd = False
    if isclose(slope, 0, abs_tol= atol):
        ctd = True
    segments[seg_num] = {'bottle': ctd, "depth": math.nan, "usable": False, "points" : []}
    return segments

def to_float(x):
    try:
        val = float(x)
    except:
        val = math.nan
    return val


ctd = '/Volumes/GeringSSD/GU201905_CTD/ctd001.evl'
with open('/Volumes/GeringSSD/eDNA_cast.txt', 'r') as cast_file:
     cast_data = cast_file.readlines()
cast_dic = {}
for cast in cast_data[1:]:
    num, depths = cast.split(maxsplit=1)
    depths = [to_float(x) for x in depths.strip().split()]
    cast_dic[int(num)] = depths
casts = cast_dic[1]


depth_data = processed_data.read_evl("", 18000, ctd)
x=depth_data.ping_time
y=depth_data.data

x_hat = [num.astype("float") for num in x]
# can adjust to get best reading possible
y_hat = savgol_filter(y, 23, 1)
rtol = 0.5
atol = 0.0002

seg_num = 0
segments = {}
current_slope = calc_slope(x_hat[0], y_hat[0], x_hat[1], y_hat[1])
segments = new_segment(segments, seg_num, current_slope)

for j in range(len(x)-1):
    new_slope = calc_slope(x_hat[j], y_hat[j], x_hat[j+1], y_hat[j+1])
    #print(new_slope)
    x_j = np.datetime_as_string(x[j])
    y_j = float(y[j])
    if not math.isclose(new_slope, current_slope, rel_tol = rtol, abs_tol=atol):
        seg_num += 1
        segments = new_segment(segments, seg_num, new_slope) # makes empty new segment
        if isclose(current_slope, 0, rel_tol = rtol, abs_tol= atol):
            segments[seg_num-1]["points"].append((x_j, y_j)) # point still in old segment
        else:
            segments[seg_num]["points"].append((x_j, y_j)) # new segment with point
        current_slope = new_slope
    else:
        segments[seg_num]["points"].append((x_j, y_j)) # puts point in current segment


combined_segs = {}
combined_segs[0] = segments[0]
curr_seg = combined_segs[0]
num_segs = 1
for num in range(1,len(segments)):
    new_seg = segments[num]
    if new_seg["bottle"] == curr_seg["bottle"]:
        curr_seg["points"] += new_seg["points"]
    else:
        combined_segs[num_segs] = new_seg
        curr_seg = combined_segs[num_segs]
        num_segs += 1

for num in range(len(combined_segs)):
    seg = combined_segs[num]
    if seg["bottle"] == True:
        times, depths = zip(*seg["points"])
        seg["depth"] = mean(depths)
        for cast_depth in casts:
            if isclose(seg["depth"], cast_depth, abs_tol=2) and cast_depth > 4.5:
                seg["usable"] = True 


for num in combined_segs:
    col = 'red'
    if combined_segs[num]["bottle"]:
        col = 'blue'
    x_seg, y_seg = zip(*combined_segs[num]["points"])
    x_seg = [np.datetime64(date) for date in x_seg]
    plt.plot(x_seg, y_seg, color = col)
plt.show()
plt.close()
with open('/Volumes/GeringSSD/segment_test.txt', 'w') as outfile:
    json.dump(combined_segs, outfile)