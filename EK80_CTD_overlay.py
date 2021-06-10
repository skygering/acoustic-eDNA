from matplotlib.pyplot import show
from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
from echolab2.processing import line
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np

def echo_plot(ek, fig, frequency, form):
    raw_data_list = ek.get_channel_data(frequencies=38000)
    raw_data = raw_data_list[frequency][0]
    for file in raw_data_list[frequency][1:]:
        raw_data.append(file)

    if form == "Pw":
        processed_power = raw_data.get_power()
        return echogram.Echogram(fig, processed_power)

    if form == "Sv":
        cal_obj = raw_data.get_calibration()
        Sv = raw_data.get_Sv(calibation=cal_obj, return_depth=True)
        return echogram.Echogram(fig, Sv, threshold=[-70,-34])


#from matplotlib.pyplot import figure, show, subplots_adjust, get_cmap
ctd_path = '/Volumes/GeringSSD/GU201905_CTD/'
# this is EK60 data recorded with EK80 software
raw_path = '/Volumes/GeringSSD/GU1905_Acoustic/EK60/'
evl_raw_file = ctd_path + "/evl_raw_matches.list"

# Read in lines of evl_raw_matches.list 
evl_raw_dic={}
with open(evl_raw_file, 'r') as file_matches:
        evl_raw = file_matches.readlines()

# Make dictionary with evl file keys and raw file values
for rows in evl_raw[23:]: # should be 1 just for testing
    evl, raw = rows.split(maxsplit=1)
    evl_raw_dic[evl] = raw.split()
    
ek80 = EK80.EK80()
for ctd in evl_raw_dic:
    # RAW FILES
    raw_infiles = []
    for raw in evl_raw_dic[ctd]:
        raw_infiles.append(raw_path + raw)

    fig, ax = plt.subplots()
    ek80.read_raw(raw_infiles)
    echogram38 = echo_plot(ek80, fig, 38000, "Sv")

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

    dt, depth = zip(*depth_data.items())
    depth_arr = np.array(depth)

    n_pings = int(evl_data[1])
    dt_arr = np.empty((n_pings), dtype='datetime64[ms]')
    for i in range(n_pings):
        dt_arr[i] = np.datetime64(dt[i])

    depth_line = line.line(ping_time = dt_arr, data = depth_arr)
    echogram38.plot_line(depth_line, linewidth=2.0)

    plt.title("Sv data with CTD depth in time order for " + ctd)
    axes = plt.gca()
    axes.set_ylim(bottom = max(depth) + 10)
    axes.set_xlim(min(dt_arr).astype('float'), max(dt_arr).astype('float'))
    echogram38.add_colorbar(fig)
    show()
    #plt.savefig(ctd_path + ctd +"Echogram.png")
    plt.close()





