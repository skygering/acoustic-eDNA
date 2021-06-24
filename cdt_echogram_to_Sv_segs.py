import find_plateau_CTD
import plot_CTD_data
import os


eDNA_cast_dic = find_plateau_CTD.read_casts('/Volumes/GeringSSD/eDNA_cast.txt')
outfile_path = '/Volumes/GeringSSD/GU201905_CTD_JSON/'

ctd_list_file = '/Volumes/GeringSSD/GU201905_CTD/CTDtoEVL.list'
evl_list = plot_CTD_data.evl_list(ctd_list_file)
for i in range(len(evl_list)):
    outfile = outfile_path + os.path.basename(evl_list[i]).replace(".evl", "") + "_segments.json"
    find_plateau_CTD.make_segments_json(eDNA_cast_dic[i+1], evl_list[i], outfile)