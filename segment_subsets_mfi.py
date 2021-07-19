'''
Load .evl/.raw file dictionary (see example in raw_overlay_ctd.py) and segment each CTD profile and take a 5 minute and 4m subset of the Sv
data around each segment where an eDNA sample was taken. Save all subsets into a nested .json file where outerlayer is cast number, the nest key is
the depth, then the frequencies, and finally the innermost layer is the Sv data points in a nested list. For each subset, calculate MFI and save image
of classification as well as a .json file with array of MFI files for each depth for each cast.

Skylar Gering - July 2021
'''

import CTD_EK_processing as process
import CTD_EK_plotting as plotting
from echolab2.instruments import EK80
import matplotlib.pyplot as plt
import os
import json
import copy

# needed paths and files
raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60/"
edna_path = "/Volumes/GeringSSD/GU1905_eDNA/"
output_path = "/Volumes/GeringSSD/GU201905_output/"
ctd_path = '/Volumes/GeringSSD/GU201905_CTD/'
ctd_list = '/CTDtoEVL.list'
evl_raw_list = output_path + "/evl_raw_matches.list"

# needed values
transducer_offset = 5

# list of needed files
asc_files = process.asc_from_list(ctd_path + ctd_list)
evl_files = process.cast_new_extension(asc_files, ".asc", ".evl")
segment_files = process.cast_new_extension(asc_files, ".asc", ".json", "segments")

# read in dictionary that matches evl files to raw files
evl_raw_dic=process.evl_raw_dic_from_file(evl_raw_list)

# read in depths of eDNA samples for each cast
eDNA_cast_dic = process.read_casts(edna_path + "eDNA_cast.txt")

ctd_subset_dic = {}
bounds_dic = {}
mfi_dic = {}

casts =  range(13,15) #range(len(eDNA_cast_dic.keys()))
for i in casts:
    cast = list(eDNA_cast_dic.keys())[i]
    mfi_depth_dic = {}
    mfi_dic[cast] = mfi_depth_dic
    file_idx = cast-1 # python counts from 0, but cast numbers start at 1y
    seg_dic = process.interactive_segment_maker(eDNA_cast_dic[cast], output_path + evl_files[file_idx], transducer_offset) # make segments
    raw_files = [os.path.normpath(raw_path + "/"+ raw) for raw in evl_raw_dic[evl_files[file_idx]]]
    subset_Sv_dic = process.interactive_subset_maker(seg_dic, raw_files, transducer_offset, toffsets = (5, 5), doffsets = (2,2)) # make subsets
    ctd_subset_dic[cast], bounds_dic[cast]  = process.subset_to_json(copy.deepcopy(subset_Sv_dic)) # save subsets in nested dictionary
    for j in range(len(subset_Sv_dic.keys())):
        sample = list(subset_Sv_dic.keys())[j]
        mfi = process.calc_MFI(subset_Sv_dic[sample])
        mfi_depth_dic[sample] = mfi.data.tolist()
        fig, ax = plt.subplots(figsize=(15,8))
        mfi_image = plotting.plot_MFI(ax, mfi, "Cast " + str(i+1) + ", Depth " + str(sample) + " : MFI Predictions")
        plt.savefig(output_path + "ctd_" +  str(i+1) + "_" + str(sample) + "_mfi.png")
        plt.close()

with open(output_path + "cast14_15_subsets.json", 'w') as outfile: # save subset dictionary to .json file
        json.dump(ctd_subset_dic, outfile)

with open(output_path + "cast14_15_subset_bounds.json", 'w') as outfile: # save subset dictionary to .json file
        json.dump(bounds_dic, outfile)

with open(output_path + "cast14_15_mfi.json", 'w') as outfile: # save subset dictionary to .json file
        json.dump(mfi_dic, outfile)





