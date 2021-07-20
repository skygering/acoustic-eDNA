import json
import CTD_EK_processing as process
import CTD_EK_plotting as plotting
import matplotlib.pyplot as plt
output_path = "/Volumes/GeringSSD/GU201905_output/"
mfi_file = "cast14_15_mfi.json"
mfi_global_file = "cast14_15_mfi_global.json"
bounds_file = "cast14_15_subset_bounds.json"
subset_file = "cast14_15_subsets.json"


subset_dic = json.load(open(output_path + subset_file))
bounds_dic = json.load(open(output_path + bounds_file))

for cast in subset_dic:
    fig, ax = plt.subplots(4, figsize = (12, 8), constrained_layout = True)
    for depth in subset_dic[cast]:
        fig.suptitle("MFI: Cast " + str(cast) + ", Depth " + str(depth))
        Sv = subset_dic[cast][depth]
        mfi_local = process.calc_MFI(Sv)
        plotting.plot_MFI(ax[0], mfi_local, "Local Norm, All Frequencies")
        mfi_local_3f = process.calc_MFI(Sv, bad_fq=[200000])
        plotting.plot_MFI(ax[1], mfi_local_3f, "Local Norm, Exclude 200kHz")
        mfi_global = process.calc_MFI(Sv, global_norm=True)
        plotting.plot_MFI(ax[2], mfi_global, "Global Norm, All frequencies")
        mfi_global_3f = process.calc_MFI(Sv, bad_fq=[200000], global_norm=True)
        plotting.plot_MFI(ax[3], mfi_global_3f, "Global Norm, Exclude 200kHz")
    plt.show()
    plt.close()