'''
Start with .asc files from CTD data and .raw files from echosounder (EK80 software). This script converts the .asc files
into .evl files and saves them at output_path. CTD traces are plotted and saved as well. Additionally, the .evl files are
matched to corresponding .raw files by time. Then, the CTD profile is overplayed on the .raw data, plotting an echogram with
the CTD profile in the foreground. Images with all 4 availible freqencies (18kHz, 38kHz, 120kHz, and 200kHz) are plotted in
a grid and saved as well. If you have more than 4 frequencies you will need to change the subplot layout.

Skylar Gering - July 2021
'''

import CTD_EK_processing as process
import CTD_EK_plotting as plotting
from echolab2.instruments import EK80
import matplotlib.pyplot as plt
import os

# plot aesthetics
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['figure.titlesize'] = 21
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

# needed paths and files
raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60/"
ctd_path = '/Volumes/GeringSSD/GU201905_CTD/'
ctd_list = '/CTDtoEVL.list'
output_path = "/Volumes/GeringSSD/GU201905_output/"

# creates a list of the .asc CTD files names
asc_files = process.asc_from_list(ctd_path + ctd_list)

# creates .evl files for each cast in CTDtoEVL.list and saves them in ctd_path
process.asc_to_evl(ctd_list, ctd_path, output_path)

# create a list of the .evl file names, which were created in the above call
evl_files = process.cast_new_extension(asc_files, ".asc", ".evl")

# gets a list of the .raw files by giving a location and then finding all the .raw files in that directory
raw_files = process.glob.glob(os.path.normpath(raw_path + "/*.raw")) # filepath in the filenames

# finds overlapping .evl and .raw files and saves the list to a dictionary (also saves to file if given outfile_path)
evl_raw_dic = process.match_raw_evl(evl_files, raw_files, evl_inpath=output_path, outfile_path=output_path)

# loop through all of the evl files
for evl in evl_raw_dic:
    # plot CTD profile
    cast_num = int(evl.replace("ctd", "").replace(".evl", ""))
    fig, ax = plt.subplots(figsize = (10,8))
    plotting.plot_evl(ax, evl, output_path, "CTD Profile: Cast " + str(cast_num))
    plt.savefig(output_path + evl.replace(".evl", ".png"))
    plt.close()
    
    # plot echograms with CTD profile overlayed
    ek80 = EK80.EK80()
    raw_infiles = [os.path.normpath(raw_path + "/"+ raw) for raw in evl_raw_dic[evl]]
    ek80.read_raw(raw_infiles)
    fig, axs = plt.subplots(2,2, figsize=(12,10), constrained_layout = True)
    fig.suptitle("Echogram with CTD Profile: Cast " + str(cast_num))
    # plot each frequency on a seperate subplot
    fq_ax = [(18000, axs[0,0]), (38000, axs[0,1]), (120000, axs[1,0]), (200000, axs[1,1])]
    for fq in fq_ax:
        echo_plot = plotting.plot_echo(fq[1], ek80, fq[0], transducer_offset= 5)
        plotting.plot_evl_trace(fq[1], echo_plot, evl, output_path)
    plt.savefig(output_path + evl.replace(".evl", "_echograms.png"))
    plt.close()




#depth_data = line.read_evl(ctd_path + evl_list[cast-1])
#x=depth_data.ping_time
#y=depth_data.data
#seg_dic = process.create_segments_dic(x, y, 0.00011)
#seg_dic = process.mark_usable_depth(seg_dic, eDNA_cast_dic[cast],5)
#plotting.plot_segments(seg_dic, "Cast 14: Segment Seperation")
#plt.show()
#plt.close()

#subset_Sv_dic= process.subset_segments_Sv(seg_dic, raw_infiles, [0, 5], [2, 2], 5)

#subset_Sv_56 = subset_Sv_dic[56] # pull out a specific depth

#mfi = process.calc_MFI(subset_Sv_56, local_norm=True)

#fig, ax = plt.subplots(figsize=(15,8))
#mfi_image = plotting.plot_MFI(ax, mfi, "Cast 14, Depth 56: MFI Predictions")
#plt.show()
#plt.close()


#for i in range(13, 14):
#    print(ctd_path + evl_list[i])
#    outfile = outfile_path + os.path.basename(json_list[i])
#    json_dic = process.interactive_segment_maker(eDNA_cast_dic[i+1], ctd_path + evl_list[i], 5)
#    print(json_dic)
