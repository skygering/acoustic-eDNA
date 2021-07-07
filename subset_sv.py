import json
import numpy as np
from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
import echolab2.processing.processed_data as pd
import os
from math import isclose
import matplotlib.pyplot as plt
import copy
import CTD_EK_processing as process

def update_Sv(Sv, xlim, ylim):
    Sv.data = Sv.data[xlim[0]:xlim[1]+1, ylim[0]:ylim[1]+1]
    Sv.ping_time = Sv.ping_time[xlim[0]:xlim[1]+1]
    Sv.depth = Sv.depth[ylim[0]:ylim[1]+1]
    Sv.transducer_draft = Sv.transducer_draft[xlim[0]:xlim[1]+1]
    Sv.n_pings = len(Sv.ping_time)
    return Sv


transducer_offset = 5
frequency_list = [18000, 38000, 120000, 200000]
toffset = 1 #needs to be integers...maybe should do seconds
doffset = 1

ctd_path = '/Volumes/GeringSSD/GU201905_CTD'
ctd_list = '/CTDtoEVL.list'
evl_raw_list = ctd_path + "/evl_raw_matches.list"
raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60"
json_path = '/Volumes/GeringSSD/GU201905_CTD_JSON/'

# get a list of asc files
asc_files = process.asc_from_list(ctd_path + ctd_list)

# This should be a function ...  it is in both overlay and here - GET RAW FILES
evl_raw_dic=process.evl_raw_dic_from_file(evl_raw_list)

# Start here
evl_files = process.cast_new_extension(asc_files, ".asc", ".evl")
json_files = process.cast_new_extension(asc_files, ".asc", ".json", "segments")

json_dic = {}
for i in range(len(json_files)):
    depth_dic = {}
    json_dic[evl_files[i].replace(".evl", "")] = depth_dic
    json_file = open(os.path.normpath(json_path + json_files[i]))
    segment_dic = json.load(json_file)

    raw_file = evl_raw_dic[evl_files[i]]
    ek80 = EK80.EK80()
    ek80.read_raw([os.path.normpath(raw_path + "/" + file) for file in raw_file])

    seg_n = 0
    for num in segment_dic:
        seg = segment_dic[num]
        if seg["usable"]:
            y_mean = float(seg["depth"])
            fq_dic = {}
            depth_dic[round(y_mean)] = fq_dic
            points_time, points_depth = zip(*seg["points"])
            points_time = [np.datetime64(x) for x in points_time]
            y_min = y_mean - doffset
            y_max = y_mean + doffset
            x_min = min(points_time) - np.timedelta64(toffset, 'm') # shoud this be from the median??
            x_max = max(points_time) + np.timedelta64(toffset, 'm')

            for frequency in frequency_list:
                raw_data_list = ek80.get_channel_data(frequencies=frequency)
                raw_data = raw_data_list[frequency][0] # I already have this in plotting code... pull out and make seperate?
                for file in raw_data_list[frequency][1:]: # could send in a list of frequencies, transducer offset, and list of files?
                    raw_data.append(file)
                cal_obj = raw_data.get_calibration()
                cal_obj.transducer_offset_z = transducer_offset
                Sv = raw_data.get_Sv(calibration=cal_obj, return_depth=True, resample_interval = raw_data.RESAMPLE_LONGEST)

                depths = np.array(Sv.depth)
                y_min_idx = np.argmin(np.abs(depths - y_min))
                y_max_idx = np.argmin(np.abs(depths - y_max))

                times = np.array(Sv.ping_time)
                x_min_idx = np.argmin(np.abs(times - x_min))
                x_max_idx = np.argmin(np.abs(times - x_max))

                subset_Sv = update_Sv(copy.deepcopy(Sv), (x_min_idx, x_max_idx), (y_min_idx, y_max_idx))
                fq_dic[str(frequency)] = subset_Sv.data.tolist()
            seg_n += 1
                
                #fig, axs = plt.subplots()
                #axs.set_title("File:" + json_files[i] + "; Segment: " + str(round(y_mean)) + "; Frequency: " + str(frequency))
                #echogram.Echogram(fig, subset_Sv, threshold=[-90,-20])

outfile_path = os.path.normpath(json_path + "/segment_subset.json")
with open(outfile_path, 'w') as outfile: # might be better to only save the  "usable" segments and then note the min/max depth and min/max time
        json.dump(json_dic, outfile)







    




