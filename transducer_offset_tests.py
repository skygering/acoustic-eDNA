from matplotlib.pyplot import show
from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
import matplotlib.pyplot as plt
import numpy as np

fig, axs = plt.subplots(2,2, figsize=(12,10))
fig.suptitle("Sv data with different transducer offset methods")

# Reading in data
raw_files = '/Volumes/GeringSSD/GU1905_Acoustic/EK60/GU19_05-D20191018-T114440.raw' # switch filepath to personal path
ek80 = EK80.EK80()
ek80.read_raw(raw_files)
raw_data_list = ek80.get_channel_data(frequencies=38000) # get data from 38000 frequency
raw_data = raw_data_list[38000][0]

# first method - original without offset attempt
cal_obj1 = raw_data.get_calibration()
Sv1 = raw_data.get_Sv(calibration=cal_obj1, return_depth=True)
print("cal 1: Original no offset")
print(cal_obj1)
echo1 = echogram.Echogram(axs[0,0], Sv1, threshold=[-80,-20])
axs[0,0].set_ylim(bottom = 25) # I have been working with this echogram so I know the depth
axs[0,0].set_title("Original - no offset attempt")

# second method - Offset attempt with transducer_offset_z
cal_obj2 = raw_data.get_calibration()
cal_obj2.transducer_offset_z = 4.5
Sv2 = raw_data.get_Sv(calibration=cal_obj2, return_depth=True)
print("Cal 2: Offset attempt with transducer_offset_z")
print(cal_obj2) # you can see that transducer_offset_z is successfully set
echo2 = echogram.Echogram(axs[0,1], Sv2, threshold=[-80,-20])
axs[0,1].set_ylim(bottom = 25)
axs[0,1].set_title("Offset attempt with transducer_offset_z")

# third method - Offset attempt with drop_keel_offset
cal_obj3 = raw_data.get_calibration()
cal_obj3.drop_keel_offset = 4.5
cal_obj3.transducer_mounting = "DropKeel"
Sv3 = raw_data.get_Sv(calibration=cal_obj3, return_depth=True)
print("Cal 3: Offset attempt with drop_keel_offset")
print(cal_obj3) # drop_keel_offset and transducer_mounting are sucessfully set
echo3 = echogram.Echogram(axs[1,0], Sv3, threshold=[-80,-20])
axs[1,0].set_ylim(bottom = 25)
axs[1,0].set_title("Offset attempt with drop_keel_offset")

# fourth method - Sv with transducer_draft set
cal_obj4 = raw_data.get_calibration()
Sv4 = raw_data.get_Sv(calibration=cal_obj4, return_depth=False)
Sv4.transducer_draft = np.full(len(Sv4.ping_time), 4.5) # sets transducer depth to 7.5 manually
Sv4.to_depth(cal_obj4) # runs the to_depth() method in processed_data.py to do offset
print("Sv4: Sv with transducer_draft set")
print(Sv4) # transducer_draft is set correctly - can only see size
print("Transducer_draft values:")
print(Sv4.transducer_draft) # this is a np array with an entry for each ping
echo4 = echogram.Echogram(axs[1,1], Sv4, threshold=[-80,-20])
axs[1,1].set_ylim(bottom = 25)
axs[1,1].set_title("Offset attempt with transducer_draft and to_depth()")

plt.show()
plt.close()