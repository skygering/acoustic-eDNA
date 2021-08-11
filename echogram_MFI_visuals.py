'''
Visuals for Skylar Gering Hollings Scholarship presentation - focus on methods for Cast 14, depth 56

Skylar Gering July 2021
'''

import CTD_EK_processing as process
import CTD_EK_plotting as plotting
import os
import matplotlib.pyplot as plt
from echolab2.instruments import EK80
from echolab2.processing import line, processed_data
from echolab2.plotting.matplotlib import echogram
import json

plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['figure.titlesize'] = 30
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

outfile_path = "/Volumes/GeringSSD/Presentation_Images/"

echo_threshold = [-75, -20]

# plots echogram with CTD overlay so I could zoom in and estimate time offset in CTD 12 and 19
raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60/"
ctd_path = '/Volumes/GeringSSD/GU201905_CTD/'
output_path = "/Volumes/GeringSSD/GU201905_output/"
ctd_list = '/CTDtoEVL.list'
ctd_num = 13 # this is cast 14 since python counts from 0
evl_raw_list = output_path + "/evl_raw_matches.list"
asc_files = process.asc_from_list(ctd_path + ctd_list)
evl_list = process.cast_new_extension(asc_files, ".asc", ".evl")
raw_list = process.glob.glob(os.path.normpath(raw_path + "/*.raw"))
evl_raw_dic=process.evl_raw_dic_from_file(evl_raw_list)


# Cast 14 Data
depth_data = line.read_evl(output_path + evl_list[13])
x=depth_data.ping_time
y=depth_data.data

# CTD Profile
fig, ax = plt.subplots(figsize=(10, 8),  constrained_layout = True)
plotting.plot_evl(ax, evl_list[13], output_path, title = "Cast 14: CTD Profile")
ax.tick_params(axis='x', labelrotation=15)
plt.savefig(outfile_path + "cast14_evl.png")
plt.close()

# 4 Frequencies Echogram
fig, ax = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
ek80 = EK80.EK80()
raw_infiles = [os.path.normpath(raw_path + "/"+ raw) for raw in evl_raw_dic[evl_list[13]]]
ek80.read_raw(raw_infiles)

fig.suptitle("Cast 14 Echograms")
plotting.plot_echo(ax[0,0], ek80, 18000, title = "18 kHz", fq_thresholds = echo_threshold, transducer_offset= 5)
ax[0,0].set_ylim(70, 0)
plotting.plot_echo(ax[0,1], ek80, 38000, title = "38 kHz", fq_thresholds = echo_threshold, transducer_offset= 5)
ax[0,1].set_ylim(70, 0)
plotting.plot_echo(ax[1,0], ek80, 120000, title = "120 kHz", fq_thresholds = echo_threshold, transducer_offset= 5)
ax[1,0].set_ylim(70, 0)
plotting.plot_echo(ax[1,1], ek80, 200000, title = "200 kHz", fq_thresholds = echo_threshold, transducer_offset= 5)
ax[1,1].set_ylim(70, 0)

plt.savefig(outfile_path + "cast14_echo.png")
plt.close()
# 38 kHz Alone
fig, ax = plt.subplots(figsize=(10, 9))
echo = plotting.plot_echo(ax, ek80, 38000, title = "Cast 14: 38 kHz", fq_thresholds = echo_threshold, transducer_offset= 5)
ax.set_ylim(70, 0)
plt.savefig(outfile_path + "cast14_38000.png")

plotting.plot_evl_trace(ax, echo, evl_list[13], output_path, time_offset=[0,0], lwidth=3)
plt.savefig(outfile_path + "cast14_38000_evl.png")
plt.close()

# segments
eDNA_cast_dic = process.read_casts('/Volumes/GeringSSD/GU1905_eDNA/eDNA_cast.txt')
json_list = process.cast_new_extension(asc_files, ".asc", ".json", "segments")
seg_dic = process.create_segments_dic(x, y, 0.00011)
seg_dic = process.mark_usable_depth(seg_dic, eDNA_cast_dic[14],5)
plotting.plot_segments(seg_dic, "Cast 14: Segment Seperation", show=False)
plt.savefig(outfile_path + "cast14_seg_seperation.png")
plt.close()


subset_Sv_dic= process.subset_segments_Sv(seg_dic, raw_infiles, [5, 5], [2, 2], 5)
Sv_18 = subset_Sv_dic["56"]["18000"]

fig, ax = plt.subplots(figsize=(15, 3),  constrained_layout = True)
echogram.Echogram(ax, Sv_18, threshold=echo_threshold)
ax.set_title("Cast 14: 18kHz Subset")
plt.savefig(outfile_path + "cast14_18000_subset.png")
plt.close()

Sv_38 = subset_Sv_dic["56"]["38000"]
fig, ax = plt.subplots(figsize=(15, 3),  constrained_layout = True)
echogram.Echogram(ax, Sv_38, threshold=echo_threshold)
ax.set_title("Cast 14: 38kHz Subset")
plt.savefig(outfile_path + "cast14_38000_subset.png")
plt.close() 

output_path = "/Volumes/GeringSSD/GU201905_output/"
mfi_file = "mfi.json"
bounds_file = "subset_bounds.json"
subset_file = "subsets.json"
mfi_dic = json.load(open(output_path + mfi_file))
subset_dic = json.load(open(output_path + subset_file))
bounds_dic = json.load(open(output_path + bounds_file))

# mfi and subset data
mfi_data = mfi_dic["14"]["56"]
bounds = bounds_dic["14"]["56"]
sub38 = subset_dic["14"]["56"]["38000"]
mfi_obj = process.processed_data_from_dic(mfi_data, bounds)

# all cast 14 depth 56 subsets
fig, ax = plt.subplots(4, figsize=(10, 8), constrained_layout=True)
echogram.Echogram(ax[0], process.processed_data_from_dic(subset_dic["14"]["56"]["18000"], bounds), threshold=echo_threshold)
ax[0].set_title("18 kHz")
echogram.Echogram(ax[1], process.processed_data_from_dic(subset_dic["14"]["56"]["38000"], bounds), threshold=echo_threshold)
ax[1].set_title("38 kHz")
echogram.Echogram(ax[2], process.processed_data_from_dic(subset_dic["14"]["56"]["120000"], bounds), threshold=echo_threshold)
ax[2].set_title("120 kHz")
echogram.Echogram(ax[3], process.processed_data_from_dic(subset_dic["14"]["56"]["200000"], bounds), threshold=echo_threshold)
ax[3].set_title("200 kHz")
fig.suptitle("Cast 14 Echograms Subset at 56m")
plt.savefig(outfile_path + "cast14_56_subsets.png")
plt.close()

# MFI for cast 14, depth 56
fig, ax = plt.subplots(figsize=(15, 3), constrained_layout = True)
plotting.plot_MFI(ax, mfi_obj, "Cast 14, Depth 56 MFI")
plt.savefig(outfile_path + "cast14_56_mfi.png")
plt.close()

fig, ax = plt.subplots(figsize=(15, 3), constrained_layout = True)
mfi_mask = process.mask_mfi(mfi_obj.data, mfi_obj.data, [[0, 0.4], [0.8, 1]])
mfi_mask_obj = process.processed_data_from_dic(mfi_mask, bounds)
plotting.plot_MFI(ax, mfi_mask_obj, "Cast 14, Depth 56 MFI Mask") # does not work
plt.savefig(outfile_path + "cast14_56_mfi_mask.png")
plt.close()

fig, ax = plt.subplots(figsize=(15, 3), constrained_layout = True)
sv_mask = process.mask_mfi(mfi_obj.data, process.processed_data_from_dic(subset_dic["14"]["56"]["38000"], bounds).data, [[0, 0.4], [0.8, 1]])
echo_plot = echogram.Echogram(ax, process.processed_data_from_dic(sv_mask, bounds), threshold=echo_threshold)
ax.tick_params(axis='x', labelrotation=15)
echo_plot.add_colorbar(fig)
ax.set_title("Mask applied to 38kHz Sv Data")
plt.savefig(outfile_path + "cast14_56_Sv_mask.png")
plt.close()

