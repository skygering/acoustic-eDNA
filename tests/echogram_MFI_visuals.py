import CTD_EK_processing as process
import CTD_EK_plotting as plotting
import os
import matplotlib.pyplot as plt
from echolab2.instruments import EK80
from echolab2.processing import line, processed_data
from echolab2.plotting.matplotlib import echogram

plt.rcParams['axes.labelsize'] = 22
plt.rcParams['axes.titlesize'] = 25
plt.rcParams['figure.titlesize'] = 28
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

outfile_path = "/Volumes/GeringSSD/Presentation_Images/"

echo_threshold = [-75, -20]

# plots echogram with CTD overlay so I could zoom in and estimate time offset in CTD 12 and 19
raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60/"
ctd_path = '/Volumes/GeringSSD/GU201905_CTD/'
ctd_list = '/CTDtoEVL.list'
ctd_num = 13 # this is cast 14 since python counts from 0
evl_raw_list = ctd_path + "/evl_raw_matches.list"
asc_files = process.asc_from_list(ctd_path + ctd_list)
evl_list = process.cast_new_extension(asc_files, ".asc", ".evl")
raw_list = process.glob.glob(os.path.normpath(raw_path + "/*.raw"))
evl_raw_dic=process.evl_raw_dic_from_file(evl_raw_list)


# Cast 14
fig, ax = plt.subplots(figsize=(10, 8),  constrained_layout = True)
plotting.plot_evl(ax, evl_list[13], ctd_path, title = "Cast 14: CTD Profile")
ax.tick_params(axis='x', labelrotation=15)
plt.savefig(outfile_path + "cast14_evl.png")
plt.close()

fig, ax = plt.subplots(figsize=(10, 8),  constrained_layout = True)
ek80 = EK80.EK80()
raw_infiles = [os.path.normpath(raw_path + "/"+ raw) for raw in evl_raw_dic[evl_list[13]]]
ek80.read_raw(raw_infiles)

# 18kHz
echo = plotting.plot_echo(ax, ek80, 18000, title = "Cast 14: 18000 kHz", 
                          fq_thresholds = echo_threshold, transducer_offset= 5)
ax.set_ylim(70, 0)
plt.savefig(outfile_path + "cast14_18000.png")

plotting.plot_evl_trace(ax, echo, evl_list[13], ctd_path)
plt.savefig(outfile_path + "cast14_18000_evl.png")
plt.close()
# 38 kHz
fig, ax = plt.subplots(figsize=(10, 8),  constrained_layout = True)
echo = plotting.plot_echo(ax, ek80, 38000, title = "Cast 14: 38000 kHz", fq_thresholds = echo_threshold, transducer_offset= 5)
ax.set_ylim(70, 0)
plt.savefig(outfile_path + "cast14_38000.png")

plotting.plot_evl_trace(ax, echo, evl_list[13], ctd_path)
plt.savefig(outfile_path + "cast14_38000_evl.png")
plt.close()


eDNA_cast_dic = process.read_casts('/Volumes/GeringSSD/eDNA_cast.txt')
json_list = process.cast_new_extension(asc_files, ".asc", ".json", "segments")
depth_data = line.read_evl(ctd_path + evl_list[13])
x=depth_data.ping_time
y=depth_data.data
seg_dic = process.create_segments_dic(x, y, 0.00011)
seg_dic = process.mark_usable_depth(seg_dic, eDNA_cast_dic[14],5)
plotting.plot_segments(seg_dic, "Cast 14: Segment Seperation", show=False)
plt.savefig(outfile_path + "cast14_seg_seperation.png")
plt.close()


subset_Sv_dic= process.subset_segments_Sv(seg_dic, raw_infiles, [0, 5], [2, 2], 5)
Sv_18 = subset_Sv_dic[56][18000]

fig, ax = plt.subplots(figsize=(10, 8),  constrained_layout = True)
echogram.Echogram(ax, Sv_18, threshold=echo_threshold)
ax.set_title("Cast 14, Frequency 18kHz Subset")
plt.savefig(outfile_path + "cast14_18000_subset.png")
plt.close()

Sv_38 = subset_Sv_dic[56][38000]
fig, ax = plt.subplots(figsize=(10, 8),  constrained_layout = True)
echogram.Echogram(ax, Sv_38, threshold=echo_threshold)
plt.savefig(outfile_path + "cast14_38000_subset.png")
plt.close() 





