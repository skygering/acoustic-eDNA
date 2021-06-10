# -*- coding: utf-8 -*-

##############################################################################
# test_EK80.py
#
# selected lines from Rick Towler's Simple_EK60_test.py
# look to his program for most of the comments
#
# jech
##############################################################################

from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from matplotlib.pyplot import figure, show, subplots_adjust, get_cmap

# this is EK60 data recorded with EK80 software
infile = '/Volumes/GeringSSD/GU1905_Acoustic/EK60/\
GU19_05-D20191018-T112417.raw'


EK_data = EK80.EK80()
EK_data.read_raw(infile)
print(EK_data)
# to get all the attributes of the EK_data (EK80) object/class you can do:
# dir(EK_data)
# vars(EK_data)
# EK_data.__dict__, or EK_data.__dict__.keys()
# the EK80 class has a frequency map, which is a dictionary keyed by frequency,
# and the value is the channel ID
freq_map = getattr(EK_data, 'frequency_map')
EK_data_freqs = list(freq_map.keys())
# build a dictinary keyed by frequency of the data we want
# if you want all the frequencies, use EK_data_freqs, otherwise make your own
# list, e.g., 
#fq_list = [18000, 38000]
#Sv_data = {k: None for k in fq_list}
#bottom_lines = {k: None for k in fq_list}
Sv_data = {k: None for k in EK_data_freqs}
bottom_lines = {k: None for k in EK_data_freqs}

# I usually label the data by frequency, rather than channel ID, but the data
# are indexed by channel ID, so need to do some translations
for fq in Sv_data.keys():
    print('Getting ', fq/1000, ' kHz data')
    tmp = EK_data.get_channel_data(frequencies = fq)
#    cal_obj = tmp[fq][0].get_calibration()
#    Sv_data[fq] = tmp[fq][0].get_Sv(cal_obj)
    Sv_data[fq] = tmp[fq][0].get_Sv()

Sv_018 = Sv_data[18000]
fig1 = plt.figure()

eg = echogram.Echogram(fig1, Sv_018, cmap='viridis')
plt.show()


# from Victoria's code
#raw_data = EK_data.get_channel_data()
#Sv_data = {18000: None, 38000: None, 120000: None, 200000: None, 70000: None}
#bottom_lines = {18000: None, 38000: None, 120000: None, 200000:None, 70000:None}
#for chan_ID in raw_data:
#    print(chan_ID, raw_data[chan_ID].frequency[0])
#    if (raw_data[chan_id].frequency[0] in Sv_data.keys()):
#        # get Sv and assign it to our dictionary
#        Sv_data[raw_data[chan_id].frequency[0]] = raw_data[chan_id].get_Sv()
#        print(Sv_data[raw_data[chan_id].frequency[0]])
#        if bottom==True:
#            bottom_lines[raw_data[chan_id].frequency[0]] = raw_data[chan_id].get_bottom()
        
#masks = {18000: None, 38000: None, 120000: None, 200000:None, 70000:None}
#for freq in Sv_data.keys():
#    if Sv_data[freq] is not None:
#        # create a mask
#        masks[freq] = mask.Mask(like=Sv_data[freq])
#
#        # create a surface exclusion line at the 10 meter RANGE.
#        surf_line = line.Line(ping_time=Sv_data[freq].ping_time, data=10)
#
#        masks[freq].apply_line(surf_line, apply_above=True)
#        # use this mask to set sample data from 0.5m above the bottom downward to NaN
#        Sv_data[freq][masks[freq]] = np.nan
           
#ax = canvas.axes
#ax.clear()
#t18 = Sv_data[18000].ping_time
#r18 = Sv_data[18000].range
#ax.invert_yaxis()
#img=ax.pcolormesh(t18, r18, self.SV_data[18000].data.T,vmin=-80, vmax=-50, cmap=cmaps().ek500)



