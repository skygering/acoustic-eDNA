import CTD_EK_processing as process
import CTD_EK_plotting as plotting
import os
import matplotlib.pyplot as plt
from echolab2.instruments import EK80

raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60/"
ctd_path = '/Volumes/GeringSSD/GU201905_CTD/'
ctd_list = '/CTDtoEVL.list'
evl_raw_list = ctd_path + "/evl_raw_matches.list"

# creates a list of the .asc CTD files names
asc_files = process.asc_from_list(ctd_path + ctd_list)

# creates .evl files for each cast in CTDtoEVL.list and saves them in ctd_path
process.asc_to_evl(ctd_list, ctd_path, "/Volumes/GeringSSD/evl_test")

# create a list of the .evl file names, which were created in the above call
evl_list = process.cast_new_extension(asc_files, ".asc", ".evl")

# gets a list of the .raw files by giving a location and then finding all the .raw files in that directory
raw_list = process.glob.glob(os.path.normpath(raw_path + "/*.raw")) # filepath in the filenames

# finds overlapping .evl and .raw files and saves the list to a dictionary (also saves to file if given outfile_path)
evl_raw_dic = process.match_raw_evl(evl_list, raw_list, evl_inpath="/Volumes/GeringSSD/GU201905_CTD", outfile_path="/Volumes/GeringSSD/evl_test")

#png_list = process.cast_new_extension(evl_list, ".evl", ".png")

#fig, ax = plt.subplots()
#plotting.plot_evl(ax, evl_list[21], evl_path=ctd_path)
#plt.close()

#evl_raw_dic=process.evl_raw_dic_from_file(evl_raw_list)

#ek80 = EK80.EK80()
#raw_infiles = []
#for raw in evl_raw_dic[evl_list[21]]:
#    raw_infiles.append(raw_path + raw)
#ek80.read_raw(raw_infiles)
#fig, axs = plt.subplots(2,2, figsize=(12,10))
#fig.suptitle("Sv data with CTD depth in time order for " + evl_list[21])

#fq_ax = [(18000, axs[0,0]), (38000, axs[0,1]), (120000, axs[1,0]), (200000, axs[1,1])]
#for fq in fq_ax:
#    echo_plot = plotting.plot_echo(ek80, fq[1], fq[0], transducer_offset= 5)
#    plotting.plot_evl_trace(echo_plot, fq[1], evl_list[21], ctd_path)

#plt.show()
#plt.close()



eDNA_cast_dic = process.read_casts('/Volumes/GeringSSD/eDNA_cast.txt')
outfile_path = '/Volumes/GeringSSD/json_test/'
json_list = process.cast_new_extension(asc_files, ".asc", ".json", "segments")
print(json_list)

#evl_list = plot_CTD_data.evl_list(ctd_path+ctd_list)
for i in range(1):
    print(ctd_path + evl_list[i])
    outfile = outfile_path + os.path.basename(json_list[i])
    json_dic = process.interactive_segment_maker(eDNA_cast_dic[i+1], ctd_path + evl_list[i], 5, outfile)
print(json_dic)