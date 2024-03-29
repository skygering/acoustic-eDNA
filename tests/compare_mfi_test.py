################################################################
# File to compare different ways of calculating MFI from Sv data
# Creates an image with MFI calculated with a local normalization
# a global normalization, and both of those with 200kHz excluded.
#
# Skylar Gering July 2021
################################################################


import json
import CTD_EK_processing as process
import CTD_EK_plotting as plotting
import matplotlib.pyplot as plt
output_path = "/Volumes/GeringSSD/GU201905_output/"
bounds_file = "cast14_15_subset_bounds.json"
subset_file = "cast14_15_subsets.json"


subset_dic = json.load(open(output_path + subset_file)) # load subsets
bounds_dic = json.load(open(output_path + bounds_file)) # load subset bounds

for cast in subset_dic: # for all casts
    for depth in subset_dic[cast]: # for all depths in a cast
        fig, ax = plt.subplots(4, figsize = (12, 8), constrained_layout = True)
        fig.suptitle("MFI: Cast " + cast + ", Depth " + depth)
        Sv = subset_dic[cast][depth] # Sv data from subset
        bounds = bounds_dic[cast][depth]
        for fq in Sv:
            Sv[fq] = process.processed_data_from_dic(Sv[fq], bounds) # create pyEcholab processed data object
        mfi_local = process.calc_MFI(Sv) # calculate MFI with local normalization with all frequencies
        plotting.plot_MFI(ax[0], mfi_local, "Local Norm, All Frequencies")
        mfi_local_3f = process.calc_MFI(Sv, bad_fq=[200000]) # calculate MFI with local normalization excluding 200kHz
        plotting.plot_MFI(ax[1], mfi_local_3f, "Local Norm, Exclude 200kHz")
        mfi_global = process.calc_MFI(Sv, global_norm=True) # calculate MFI with global normalization with all frequencies
        plotting.plot_MFI(ax[2], mfi_global, "Global Norm, All frequencies")
        mfi_global_3f = process.calc_MFI(Sv, bad_fq=[200000], global_norm=True) # calculate MFI with global normalization excluding 200kHz
        plotting.plot_MFI(ax[3], mfi_global_3f, "Global Norm, Exclude 200kHz")
        plt.savefig(output_path + "compare_mfi_" + cast + depth + ".png")
        plt.close()