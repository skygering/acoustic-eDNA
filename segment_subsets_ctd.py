import CTD_EK_processing as process
import CTD_EK_plotting as plotting
from echolab2.instruments import EK80
import matplotlib.pyplot as plt
import os
import json

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

for i in  range(0,1):#range(len(eDNA_cast_dic.keys())):
    cast = list(eDNA_cast_dic.keys())[i]
    file_idx = cast-1 # python counts from 0, but cast numbers start at 1y
    seg_dic = process.interactive_segment_maker(eDNA_cast_dic[cast], output_path + evl_files[file_idx], transducer_offset)
    raw_files = [os.path.normpath(raw_path + "/"+ raw) for raw in evl_raw_dic[evl_files[file_idx]]]
    subset_Sv_dic = process.interactive_subset_maker(seg_dic, raw_files, transducer_offset, toffsets = (0, 5), 
                                                     doffsets = (2,2), outfile_path=output_path + "test.json")







