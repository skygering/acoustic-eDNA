import numpy as np
from matplotlib.pyplot import figure, show, subplots_adjust, get_cmap
from echolab2.processing import processed_data, line
from echolab2.plotting.matplotlib import echogram
from echolab2.instruments import EK60
#from echolab2.instruments import EK80
# EK60      '/Volumes/GeringSSD/D20090916-T135430.raw'
# EK80      '/Volumes/GeringSSD/GU1905_Acoustic/EK60/GU19_05-D20191015-T112548.raw'
rawfiles = ['/Volumes/GeringSSD/D20090916-T135430.raw']
#rawfiles = ['/Volumes/GeringSSD/GU1905_Acoustic/EK60/\
#GU19_05-D20191101-T173021.raw']
fig = figure()
# Set some properties for the sub plot layout.
subplots_adjust(left=0.11, bottom=0.1, right=0.98, top=.93, wspace=None,
                hspace=0.9)

ek60 = EK60.EK60()
#ek60 = EK80.EK80()
ek60.read_raw(rawfiles)
print(ek60)

raw_data_38 = ek60.get_channel_data(frequencies=38000)
raw_38 = raw_data_38[38000][0]

ax_1 = fig.add_subplot(3, 1, 1)
echogram_1 = echogram.Echogram(ax_1, raw_38, 'power')

# Playing around
echogram_1.add_colorbar(fig)

ax_1.set_title("Power as stored in raw_data object")

processed_power = raw_38.get_power()
print(processed_power)

ax_2 = fig.add_subplot(3, 1, 2)
echogram_2 = echogram.Echogram(ax_2, processed_power)
ax_2.set_title("Power data in time order")

cal_obj = raw_38.get_calibration()
Sv = raw_38.get_Sv(calibation=cal_obj, return_depth=True)
print(Sv)

ax_3 = fig.add_subplot(3, 1, 3)
echogram_3 = echogram.Echogram(ax_3, Sv, threshold=[-70,-34])
ax_3.set_title("Sv data in time order")
show()

# Create another matplotlib figure.
fig = figure()
# Set some properties for the sub plot layout.
subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=.93, wspace=None,
                hspace=0.5)

angle_cmap = get_cmap('plasma')

# Now request angles data in time order.
angles_along, angles_athwart = raw_38.get_physical_angles()
print(angles_along)
print(angles_athwart)

# Create another axis.
ax_1 = fig.add_subplot(2, 1, 1)
# Create an echogram which will display on our newly created axis.
echogram_3 = echogram.Echogram(ax_1, angles_along, cmap=angle_cmap)
ax_1.set_title("angles_alongship data in time order")

# Create another axis.
ax_2 = fig.add_subplot(2, 1, 2)
# Create an echogram which will display on our newly created axis.
echogram_3 = echogram.Echogram(ax_2, angles_athwart, cmap=angle_cmap)
ax_2.set_title("angles_athwartship data in time order")

# Show our figure.
show()
