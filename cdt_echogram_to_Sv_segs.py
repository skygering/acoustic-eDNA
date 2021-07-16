import CTD_EK_processing as process
import CTD_EK_plotting as plotting
import os
import matplotlib.pyplot as plt
from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
from echolab2.processing import line


raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60/"
ctd_path = '/Volumes/GeringSSD/GU201905_CTD/'
ctd_list = '/CTDtoEVL.list'
evl_raw_list = ctd_path + "/evl_raw_matches.list"
json_path = "/Volumes/GeringSSD/GU201905_CTD_JSON/"

cast = 14

# creates a list of the .asc CTD files names
asc_files = process.asc_from_list(ctd_path + ctd_list)

# creates .evl files for each cast in CTDtoEVL.list and saves them in ctd_path
#process.asc_to_evl(ctd_list, ctd_path, ctd_path)

# create a list of the .evl file names, which were created in the above call
evl_list = process.cast_new_extension(asc_files, ".asc", ".evl")

# gets a list of the .raw files by giving a location and then finding all the .raw files in that directory
raw_list = process.glob.glob(os.path.normpath(raw_path + "/*.raw")) # filepath in the filenames

# finds overlapping .evl and .raw files and saves the list to a dictionary (also saves to file if given outfile_path)
evl_raw_dic = process.match_raw_evl(evl_list, raw_list, evl_inpath=ctd_path, outfile_path=ctd_path)
#evl_raw_dic=process.evl_raw_dic_from_file(evl_raw_list)

ek80 = EK80.EK80()
raw_infiles = [os.path.normpath(raw_path + "/"+ raw) for raw in evl_raw_dic[evl_list[cast-1]]]
ek80.read_raw(raw_infiles)

#fig, axs = plt.subplots(2,2, figsize=(12,10))
#fig.suptitle("Sv data with CTD depth in time order for " + evl_list[cast-1])

#fq_ax = [(18000, axs[0,0]), (38000, axs[0,1]), (120000, axs[1,0]), (200000, axs[1,1])]
#for fq in fq_ax:
#    echo_plot = plotting.plot_echo(fq[1], ek80, fq[0], transducer_offset= 5)
#    plotting.plot_evl_trace(fq[1], echo_plot, evl_list[cast-1], ctd_path)

#plt.show()
#plt.close()



eDNA_cast_dic = process.read_casts('/Volumes/GeringSSD/eDNA_cast.txt') # potentially get rid of header ...
outfile_path = '/Volumes/GeringSSD/json_test/'
json_list = process.cast_new_extension(asc_files, ".asc", ".json", "segments")

#seg_dic = process.interactive_segment_maker(eDNA_cast_dic[cast], ctd_path + evl_list[cast-1], 5, outfile_path = json_path + json_list[cast-1])
#plt.close()

depth_data = line.read_evl(ctd_path + evl_list[cast-1])
x=depth_data.ping_time
y=depth_data.data
seg_dic = process.create_segments_dic(x, y, 0.00011)
seg_dic = process.mark_usable_depth(seg_dic, eDNA_cast_dic[cast],5)
#plotting.plot_segments(seg_dic, "Cast 14: Segment Seperation")
#plt.show()
#plt.close()

subset_Sv_dic= process.subset_segments_Sv(seg_dic, raw_infiles, [0, 5], [2, 2], 5)

subset_Sv_56 = subset_Sv_dic[56] # pull out a specific depth

mfi = process.calc_MFI(subset_Sv_56, local_norm=True)

fig, ax = plt.subplots()
mfi_image = plotting.plot_MFI(ax, fig, mfi)
fig.colorbar(mfi_image)
plt.show()
plt.close()


#for i in range(13, 14):
#    print(ctd_path + evl_list[i])
#    outfile = outfile_path + os.path.basename(json_list[i])
#    json_dic = process.interactive_segment_maker(eDNA_cast_dic[i+1], ctd_path + evl_list[i], 5)
#    print(json_dic)
