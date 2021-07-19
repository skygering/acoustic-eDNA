import json
import numpy as np
from echolab2.processing import processed_data
import math
from echolab2.plotting.matplotlib import echogram
import matplotlib.pyplot as plt

output_path = "/Volumes/GeringSSD/GU201905_output/"
mfi_file = "cast14_15_mfi.json"
bounds_file = "cast14_15_subset_bounds.json"
subset_file = "cast14_15_subsets.json"

mfi_dic = json.load(open(output_path + mfi_file))
subset_dic = json.load(open(output_path + subset_file))
bounds_dic = json.load(open(output_path + bounds_file))

mfi_14_56 = np.array(mfi_dic["14"]["56"])
mask = np.where((0 < mfi_14_56) & (mfi_14_56< 0.4), 1, 0)

sub_14_56_38 = np.array(subset_dic["14"]["56"]["38000"])

mask_Sv = np.multiply(mask, sub_14_56_38)
print(mask_Sv)

mask_obj = processed_data.processed_data("None", math.nan, "Sv")
mask_obj.data = mask_Sv
pings = bounds_dic["14"]["56"]["ping_time"]
mask_obj.ping_time = np.array([np.datetime64(x) for x in pings])
mask_obj.n_pings = len(pings)
mask_obj.depth = np.array(bounds_dic["14"]["56"]["depth"])

fig, ax = plt.subplots()
echogram.Echogram(ax, mask_obj, threshold=[-75, -20])
plt.show()







