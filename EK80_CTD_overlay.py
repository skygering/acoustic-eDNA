####################################################################################
# EK80_CTD_overlay.py
#
# Overlaying depth data from CTD over select echograms
# Uses Rick Towler's pyEcholab package - see files for more extensive documentation
# 
# This script assumes that you have access to both CTD data and .raw files from
# the same cruise (and specify the paths to that data below) 
# and have run match_acoustic_CTD.py and thus have evl_raw_matches.list
# 
# Skylar Gering 06/10/21
####################################################################################

from matplotlib.pyplot import show
from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
from echolab2.processing import line, processed_data
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
from echolab2.instruments.util.simrad_calibration import calibration

def echo_plot(ek, fig, frequency, form, echo_bottom, x_lims, lower_threshold = -90):
    """
    Creates an echogram with a depth profile from CTD data overlayed over time
    Inputs: EK object (EK60 or EK80 - only EK80 tested), matplotlib figure
            integer frequency of echosounder, and string "Sv" or "Pw" for
            if you want Sv data or Power data plotted
    Outputs: An echogram object
    """
    fig.set_title(str(frequency) + "Hz")
    fig.set_ylim(echo_bottom, 0)
    fig.set_xlim(x_lims[0], x_lims[1])
    raw_data_list = ek.get_channel_data(frequencies=frequency) # get data from specific frequency
    raw_data = raw_data_list[frequency][0]

    for file in raw_data_list[frequency][1:]: 
        # if there is more than one .raw file, this appends all data of specified frequency
        raw_data.append(file)

    if form == "Pw": # plot power data
        processed_power = raw_data.get_power() # get power data
        return echogram.Echogram(fig, processed_power)

    if form == "Sv": # plot Sv data
        cal_obj = raw_data.get_calibration()
        Sv = raw_data.get_Sv(calibration=cal_obj, return_depth=False) # get Sv data
        Sv.transducer_draft = np.full(len(Sv.ping_time), 4.5) # sets transducer depth to 7.5 manually - might be wrong
        Sv.to_depth(cal_obj)
        return echogram.Echogram(fig, Sv, threshold=[lower_threshold,-20])

# paths and files needed
ctd_path = '/Volumes/GeringSSD/GU201905_CTD/'
raw_path = '/Volumes/GeringSSD/GU1905_Acoustic/EK60/'
evl_raw_file = ctd_path + "/evl_raw_matches.list"
# lower frequency threshold for noise
lower_threshold = -70

# Read in lines of evl_raw_matches.list 
evl_raw_dic={}
with open(evl_raw_file, 'r') as file_matches:
        evl_raw = file_matches.readlines()

# Make dictionary with evl file keys and raw file values
for rows in evl_raw[3:5]: # should be 1 just for testing
    evl, raw = rows.split(maxsplit=1)
    evl_raw_dic[evl] = raw.split()
ek80 = EK80.EK80()

for ctd in evl_raw_dic:
    # RAW FILES
    raw_infiles = []
    for raw in evl_raw_dic[ctd]:
        raw_infiles.append(raw_path + raw)
    ek80.read_raw(raw_infiles)

    #EVL FILES
    with open(ctd_path + ctd, 'r') as evl_infile:
        evl_data = evl_infile.readlines()
    depth_data = {}
    for rows in evl_data[2:]: # exclude the first two lines which are headers
        rows.rstrip()
        # Pull out the date, time, and depth (one additional variable ignored) -  all of these are strings
        (d, t, depth, tmpVars) = rows.split(maxsplit=3)
        dt = datetime.strptime(d + t, "%Y%m%d%H%M%S%f") # use date and time strings to make datetime object
        depth_data[dt] = float(depth) # negate depth and make into float

    dt,depth= zip(*[(dt, depth) for dt, depth in depth_data.items() if depth > 0])
    depth_arr = np.array(depth)
    n_pings = len(depth)
    dt_arr = np.empty((n_pings), dtype='datetime64[ms]')
    for i in range(n_pings):
        dt_arr[i] = np.datetime64(dt[i])


    # PLOTTING
    fig, axs = plt.subplots(2,2, figsize=(10,10))
    fig.suptitle("Sv data with CTD depth in time order for " + ctd)
    echo_bottom = max(depth) + 10
    x_lim = (min(dt_arr).astype('float'), max(dt_arr).astype('float'))

    fq_ax = [(18000, axs[0,0]), (38000, axs[0,1]), (120000, axs[1,0]), (200000, axs[1,1])]
    # Plot 4x4 grid of each frequency for a CTD/echogram pairing
    for fq in fq_ax:
        echo_fq = echo_plot(ek80, fq[1], fq[0], "Sv", echo_bottom, x_lim, lower_threshold)
        depth_line = line.line(ping_time = dt_arr, data = depth_arr)
        echo_fq.plot_line(depth_line, linewidth=2.5, color = "black")
    #plt.savefig(ctd_path + ctd.replace(".evl", "") + "_echogram_" + str(lower_threshold) + ".png")
    plt.show()
    plt.close()

