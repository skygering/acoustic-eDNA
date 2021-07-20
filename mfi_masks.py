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

c = "15"
d = "9"

mfi_dic = json.load(open(output_path + mfi_file))
#mfi_3fq_dic = json.load(open(output_path + mfi_3fq_file))
subset_dic = json.load(open(output_path + subset_file))
bounds_dic = json.load(open(output_path + bounds_file))

sub_14_56_38 = subset_dic[c][d]["38000"]

mask_Sv = process.mask_mfi(mfi_dic[c][d], sub_14_56_38, [0, 0.4])
#mask_Sv_3fq = process.mask_mfi(mfi_3fq_dic[c][d], sub_14_56_38, [0, 0.4])

mask_obj = process.processed_data_from_dic(mask_Sv, bounds_dic[c][d])
#mask_3fq_obj = process.processed_data_from_dic(mask_Sv_3fq, bounds_dic[c][d])
Sv18_obj =process.processed_data_from_dic(sub_14_56_38, bounds_dic[c][d])

fig, ax = plt.subplots(2, figsize = (11, 4))
echogram.Echogram(ax[0], mask_obj, threshold=[-90, -20])
echogram.Echogram(ax[1], Sv18_obj, threshold=[-90, -20])
plt.show()

print(process.calc_ABC(mask_obj))
#print(process.calc_ABC(mask_3fq_obj))






