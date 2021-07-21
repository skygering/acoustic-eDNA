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
mfi_file = "cast14_15_mfi.json"
#mfi_3fq_file = "cast14_15_mfi_3fq.json"
bounds_file = "cast14_15_subset_bounds.json"
subset_file = "cast14_15_subsets.json"

mfi_dic = json.load(open(output_path + mfi_file))
subset_dic = json.load(open(output_path + subset_file))
bounds_dic = json.load(open(output_path + bounds_file))

for cast in mfi_dic:
    for depth in mfi_dic[cast]:
        mfi = mfi_dic[cast][depth]
        bounds = bounds_dic[cast][depth]
        sub38 = subset_dic[cast][depth]["38000"] # using 38kHz for ABC calculation
        masked_Sv = process.mask_mfi(mfi, sub38, [[0, 0.4], [0.8, 1]]) # mask for all fish (with and without swimbladder)
        # create processed data objects for ease of plotting
        mask_Sv_obj = process.processed_data_from_dic(masked_Sv, bounds)
        Sv18_obj =process.processed_data_from_dic(sub38, bounds)
        # plot comparison of Sv data and masked Sv data
        fig, ax = plt.subplots(2, figsize = (11, 4))
        echogram.Echogram(ax[0], mask_Sv_obj, threshold=[-90, -20])
        echogram.Echogram(ax[1], Sv18_obj, threshold=[-90, -20])
        plt.show()
        plt.close()

        print(process.calc_ABC(mask_Sv_obj)) # print and calculate ABC of masked Sv






