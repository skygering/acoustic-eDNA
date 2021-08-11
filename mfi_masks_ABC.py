'''
Load MFI, subset, and bounds files to create a mask for fish using MFI calcualtions.
This mask is then applied to the 38kHz Sv data to calculate the ABC value.

Skylar Gering July 2021
'''

import json
import numpy as np
from echolab2.processing import processed_data
import math
from echolab2.plotting.matplotlib import echogram
import matplotlib.pyplot as plt
import CTD_EK_processing as process

output_path = "/Volumes/GeringSSD/GU201905_output/"
mfi_file = "mfi.json"
bounds_file = "subset_bounds.json"
subset_file = "subsets.json"

mfi_dic = json.load(open(output_path + mfi_file))
subset_dic = json.load(open(output_path + subset_file))
bounds_dic = json.load(open(output_path + bounds_file))

ABC_dic = {}
for cast in mfi_dic:
    abc_depth_dic = {}
    ABC_dic[cast] = abc_depth_dic
    for depth in mfi_dic[cast]:
        mfi = mfi_dic[cast][depth]
        if mfi is not None:
            bounds = bounds_dic[cast][depth]
            sub38 = subset_dic[cast][depth]["38000"] # using 38kHz for ABC calculation
            masked_Sv = process.mask_mfi(mfi, sub38, [[0, 0.4], [0.8, 1]]) # mask for all fish (with and without swimbladder)
            # create processed data objects for ease of plotting
            mask_Sv_obj = process.processed_data_from_dic(masked_Sv, bounds)
            Sv18_obj =process.processed_data_from_dic(sub38, bounds)
            # plot comparison of Sv data and masked Sv data
            fig, ax = plt.subplots(2, figsize = (11, 4), constrained_layout = True)
            fig.suptitle("Cast " + cast + ", Depth " + depth)
            echogram.Echogram(ax[0], Sv18_obj, threshold=[-90, -20])
            ax[0].set_title("Unmasked 38kHz Sv Data")
            echogram.Echogram(ax[1], mask_Sv_obj, threshold=[-90, -20])
            ax[1].set_title("Masked 38kHz Sv Data")
            plt.show()
            plt.close()

            abc_val = process.calc_ABC(mask_Sv_obj) # print and calculate ABC of masked Sv
            abc_depth_dic[depth] = abc_val

with open(output_path + "abc.json", 'w') as outfile: # save subset dictionary to .json file
        json.dump(ABC_dic, outfile)






