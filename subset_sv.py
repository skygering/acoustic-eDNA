import json
import numpy as np
from numpy.core.fromnumeric import ndim, shape, size
from numpy.lib.function_base import append
from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
import echolab2.processing.processed_data as pd
import os
from math import isclose
import matplotlib.pyplot as plt
import copy
import CTD_EK_processing as process
import CTD_EK_plotting as plotting

def raw_to_Sv(ek, fq, transducer_offset):
    raw_data_list = ek.get_channel_data(frequencies=fq)
    raw_data = raw_data_list[fq][0] # I already have this in plotting code... pull out and make seperate?
    for file in raw_data_list[fq][1:]: # could send in a list of frequencies, transducer offset, and list of files?
        raw_data.append(file)
    cal_obj = raw_data.get_calibration()
    cal_obj.transducer_offset_z = transducer_offset
    return raw_data.get_Sv(calibration=cal_obj, return_depth=True), cal_obj


def crop_Sv(Sv, xlim, ylim):
    # data has dims (ping_time, depth) - rows are ping_time, columns are depth
    Sv.data = Sv.data[xlim[0]:xlim[1]+1, ylim[0]:ylim[1]+1] 
    Sv.ping_time = Sv.ping_time[xlim[0]:xlim[1]+1]
    Sv.depth = Sv.depth[ylim[0]:ylim[1]+1]
    Sv.transducer_offset = Sv.transducer_offset[xlim[0]:xlim[1]+1]
    Sv.n_pings = len(Sv.ping_time)
    Sv.heave = Sv.heave[xlim[0]:xlim[1]+1]
    return Sv


transducer_offset = 5


ctd_path = '/Volumes/GeringSSD/GU201905_CTD'
ctd_list = '/CTDtoEVL.list'
evl_raw_list = ctd_path + "/evl_raw_matches.list"
raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60"
json_path = '/Volumes/GeringSSD/GU201905_CTD_JSON/'
subset_path = '/Volumes/GeringSSD/Subset_Sv/'

# get a list of asc files
asc_files = process.asc_from_list(ctd_path + ctd_list)

# This should be a function ...  it is in both overlay and here - GET RAW FILES
evl_raw_dic=process.evl_raw_dic_from_file(evl_raw_list)

# Start here
evl_files = process.cast_new_extension(asc_files, ".asc", ".evl")
json_files = process.cast_new_extension(asc_files, ".asc", ".json", "segments")

json_dic = {}
for i in range(13, 15):#range(len(json_files)):
    depth_dic = {}
    ctd_name = evl_files[i].replace(".evl", "")
    json_dic[ctd_name] = depth_dic
    json_file = open(os.path.normpath(json_path + json_files[i]))
    segment_dic = json.load(json_file)

    raw_file = evl_raw_dic[evl_files[i]]
    ek80 = EK80.EK80()
    ek80.read_raw([os.path.normpath(raw_path + "/" + file) for file in raw_file])

    for num in segment_dic:
        pre_toffset = 0
        post_toffset = 5 #needs to be integers number of minutes
        top_doffset = 2
        bot_doffset = 2
        seg = segment_dic[num]
        if seg["usable"]:
            y_mean = float(seg["depth"])
            fq_dic = {}
            depth_dic[round(y_mean)] = fq_dic
            points_time, points_depth = zip(*seg["points"])
            points_time = [np.datetime64(x) for x in points_time]

            if y_mean < 11: # this is going to need more finess 
                top_doffset = 0.5
                bot_doffset = 3.5
            y_min = y_mean - top_doffset
            y_max = y_mean + bot_doffset
            x_min = min(points_time) - np.timedelta64(pre_toffset, 'm') 
            x_max = max(points_time) + np.timedelta64(post_toffset, 'm')

            Sv18, cal18 = raw_to_Sv(ek80, 18000, transducer_offset)

            # once we resample, these should be the same for each frequency
            depths = np.array(Sv18.depth)
            y_min_idx = np.argmin(np.abs(depths - y_min)) # find the index of the depth value that is closest to y_min
            y_max_idx = np.argmin(np.abs(depths - y_max))

            times = np.array(Sv18.ping_time)
            x_min_idx = np.argmin(np.abs(times - x_min))
            x_max_idx = np.argmin(np.abs(times - x_max))

            fig, axs = plt.subplots(2,2, figsize=(12,10))
            fig.suptitle("File:" + json_files[i] + "; Depth: " + str(round(y_mean)))
            fq_ax = [(18000, axs[0,0]), (38000, axs[0,1]), (120000, axs[1,0]), (200000, axs[1,1])]   

            bad_fq = []
            if y_mean > 60:
                bad_fq = append(bad_fq, 200000)
                if y_mean > 260:
                    bad_fq = append(bad_fq, 120000)

            for fq in fq_ax:
                if fq[0] not in bad_fq:
                    Sv, cal =  raw_to_Sv(ek80, fq[0], transducer_offset)
                    if cal.sample_interval != cal18.sample_interval:
                        Sv.match_samples(Sv18)
                    subset_Sv = crop_Sv(copy.deepcopy(Sv), (x_min_idx, x_max_idx), (y_min_idx, y_max_idx))
                    fq_dic[str(fq[0])] = subset_Sv.data.tolist()
                    
                    fq[1].set_title(str(fq[0]) + "Hz")
                    echogram.Echogram(fq[1], subset_Sv, threshold=[-90,-20])
            #plt.show()
            print(subset_path + ctd_name + "_" + str(round(y_mean)) + '_subset.png')
            plt.savefig(subset_path + ctd_name + "_" + str(round(y_mean)) + '_subset.png')
            plt.close()

outfile_path = os.path.normpath(subset_path + "/segment_subset_CTD14-15_2.json")
with open(outfile_path, 'w') as outfile: # might be better to only save the  "usable" segments and then note the min/max depth and min/max time
    json.dump(json_dic, outfile)







    




