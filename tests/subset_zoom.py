import json
from echolab2.instruments import EK80
import CTD_EK_processing as process
import os
import matplotlib.pyplot as plt
import numpy as np
from echolab2.plotting.matplotlib import echogram

raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60/"
ctd_path = '/Volumes/GeringSSD/GU201905_CTD/'
output_path = "/Volumes/GeringSSD/GU201905_output/"
ctd_list = '/CTDtoEVL.list'
mfi_file = "mfi.json"
bounds_file = "subset_bounds.json"
subset_file = "subsets.json"
evl_raw_list = output_path + "/evl_raw_matches.list"

mfi_dic = json.load(open(output_path + mfi_file))
subset_dic = json.load(open(output_path + subset_file))
bounds_dic = json.load(open(output_path + bounds_file))

evl_raw_dic = process.evl_raw_dic_from_file(evl_raw_list)
asc_files = process.asc_from_list(ctd_path + ctd_list)
evl_list = process.cast_new_extension(asc_files, ".asc", ".evl")
raw_list = process.glob.glob(os.path.normpath(raw_path + "/*.raw"))

mfi = mfi_dic["15"]["42"]
subset = subset_dic["15"]["42"]["38000"]
bounds = bounds_dic["15"]["42"]

subset_obj = process.processed_data_from_dic(subset, bounds)

masked_Sv = process.mask_mfi(mfi, subset, [[0, 0.4], [0.8, 1]])
mask_Sv_obj = process.processed_data_from_dic(masked_Sv, bounds)
abc_val = process.calc_ABC(mask_Sv_obj)

print(abc_val)

print(np.max(subset))
fig, ax = plt.subplots()
echogram.Echogram(ax, subset_obj, threshold=[-75, -20])
plt.show()
plt.close()

fig, ax = plt.subplots()
ax.hist(subset, bins = 20)
plt.show()
plt.close()









