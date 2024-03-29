######################################################################
# Script to test the pyEcholab line.read_evl method (written by Gering,
# modified by Towler July 2021) 
#
# Skylar Gering July 2021
######################################################################

from echolab2.processing import line
from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
import matplotlib.pyplot as plt

raw_files = '/Volumes/GeringSSD/GU1905_Acoustic/EK60/GU19_05-D20191018-T114440.raw' # switch filepath to personal path
ek80 = EK80.EK80()
ek80.read_raw(raw_files)
frequency = 38000
raw_data_list = ek80.get_channel_data(frequencies=frequency) # get data from 38000 frequency
raw_data = raw_data_list[38000][0]
# CTD data read in and made into a line
depth_line = line.read_evl("/Volumes/GeringSSD/GU201905_CTD/ctd001.evl")
# Print out looks good
print(depth_line)

# get Sv data
cal_obj = raw_data.get_calibration()
Sv = raw_data.get_Sv(calibation=cal_obj, return_depth=False)

# plotting
fig, ax = plt.subplots()
plt.ylim(bottom = 25)
plt.xlim((min(depth_line.ping_time).astype('float')-10, max(depth_line.ping_time).astype('float')))
echo_plot = echogram.Echogram(fig, Sv, threshold=[-80,-20])
echo_plot.plot_line(depth_line, linewidth=2.5, color = "black")
# This is going to show a depth trace that goes off of the edge of the echogram for the ctd001 data-
# this is because the CTD data was taken over 2-3 echograms. I only used one for testing purposes
plt.show()