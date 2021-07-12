
"""
This script test the processed_data.regird method which resamples the vertical
axis of a processed_data object to match the sample interval of the provided
processed_data object:

    # regrid the 38 kHz 512us pulse length data to the 18 kHz 1024us grid
    Sv_38_regrid.regrid(Sv_18)

so it is fairly slow for large arrays. It uses np.matmult which is about 15-20%
faster in optimized versions of Numpy (with MKL+BLAS) vs standard Numpy.

"""

import sys
import time
import numpy as np
from echolab2.instruments import EK80
import echolab2.processing.processed_data as processed_data
from echolab2.plotting.matplotlib import echogram
import matplotlib.pyplot as plt


def read_write_callback(filename, cumulative_pct, cumulative_bytes, userref):
    '''
    read_write_callback is a simple example of using the progress_callback
    functionality of the EK80.read_raw and EK80.write_raw methods.
    '''

    if cumulative_pct > 100:
        return

    if cumulative_pct == 0:
        sys.stdout.write(filename)

    if cumulative_pct % 4:
        sys.stdout.write('.')

    if cumulative_pct == 100:
        sys.stdout.write('  done!\n')



# specify the ping to plot
plot_ping = 24

#  source raw file
rawfiles = ['C:/EK Test Data/EK80/CW/bad offsets/leg 1 EK80/DY2104-D20210607-T030806.raw']
# Echoview export of original 38 kHz data
ev_export = 'C:/EK Test Data/EK80/CW/bad offsets/leg 1 EK80/DY2104-D20210607-T030806_38kHz.mat'
# Echoview export of 38 kHz data exported from Match Geometry virtual to 18 kHz 1024us grid
ev_match_geo_export = 'C:/EK Test Data/EK80/CW/bad offsets/leg 1 EK80/DY2104-D20210607-T030806_38kHz_resample_to_1024us.mat'


# Read in the EV regridded 38 (match geometry operator export)
print('Reading the echoview file %s' % (ev_match_geo_export))
ev_Sv_38_regrid = processed_data.read_ev_mat('', 38000, ev_match_geo_export, data_type='Sv')

# Read in the EV 38
print('Reading the echoview file %s' % (ev_match_geo_export))
ev_Sv_38 = processed_data.read_ev_mat('', 38000, ev_export, data_type='Sv')

# Read the raw data
ek80 = EK80.EK80()
print('Reading raw files...')
ek80.read_raw(rawfiles, frequencies=[18000,38000], progress_callback=read_write_callback)

# Get the 18 and 38 kHz raw data.
print('Getting Sv...')
raw_data = ek80.get_channel_data(frequencies=[18000,38000])
raw_data_38_0512 = raw_data[38000][0]
raw_data_18_1024 = raw_data[18000][0]

# Convert to Sv
cal_38 = raw_data_38_0512.get_calibration()
cal_18 = raw_data_18_1024.get_calibration()
Sv_38 = raw_data_38_0512.get_Sv(calibration=cal_38)
Sv_18 = raw_data_18_1024.get_Sv(calibration=cal_18)

# Proove that you can't subtract processed data objects with different axes
try:
    Sv18_minus_Sv38 = Sv_18 - Sv_38
except Exception as e:
    print("Trying to subtract processed_data objects with different vertical axes " +
           "results in the following error:")
    print("    " + str(e))


# Regrid 38 to 18 kHz vertical grid
print('Regridding...')
Sv_38_regrid = Sv_38.copy()
s_time = time.perf_counter()
Sv_38_regrid.match_samples(Sv_18)
print('Regrid elapsed time:' + str(time.perf_counter()-s_time))



# compute mean sample values to compare before/after regrid
print()
print('Mean sample values for ping ' + str(plot_ping))
Sv_38.to_linear()
val = (10.0 * np.log10(np.nansum(Sv_38.data[plot_ping,:]) / Sv_38.data.shape[1]))
print('Echolab Sv - original: ' + str(val))
Sv_38_regrid.to_linear()
val = (10.0 * np.log10(np.nansum(Sv_38_regrid.data[plot_ping,:]) / Sv_38_regrid.data.shape[1]))
print('Echolab Sv - regridded: ' + str(val))
print()
ev_Sv_38.to_linear()
val = (10.0 * np.log10(np.nansum(ev_Sv_38.data[plot_ping,:]) / ev_Sv_38.data.shape[1]))
print('Echoview Sv - original: ' + str(val))
ev_Sv_38_regrid.to_linear()
val = (10.0 * np.log10(np.nansum(ev_Sv_38_regrid.data[plot_ping,:]) / ev_Sv_38_regrid.data.shape[1]))
print('Echoview Sv - regridded: ' + str(val))

# convert back to log form for plotting
Sv_38_regrid.to_log()
Sv_38.to_log()
ev_Sv_38_regrid.to_log()
ev_Sv_38.to_log()

# Plot this all up
fig = plt.figure()
eg = echogram.Echogram(fig, Sv_38, threshold=[-70,-34])
eg.add_colorbar(fig)
eg.axes.set_title("Original - Sv 38 kHz 512us")
fig = plt.figure()
eg = echogram.Echogram(fig, Sv_38_regrid, threshold=[-70,-34])
eg.add_colorbar(fig)
eg.axes.set_title("Regridded Sv 38 kHz to match 18 kHz 1024us")

# Plot the difference between our 18 kHz at 1024us and the regridded 38 kHz data

Sv18_minus_Sv38 = Sv_18 - Sv_38_regrid
fig = plt.figure()
eg = echogram.Echogram(fig, Sv18_minus_Sv38, threshold=[-30,30])
eg.add_colorbar(fig)
eg.axes.set_title("18 kHz minus regridded 38 kHz")

fig2 = plt.figure()
plt.plot(ev_Sv_38_regrid[plot_ping], ev_Sv_38_regrid.range, label='Echoview Regrid', color='blue', linewidth=1)
plt.plot(Sv_38_regrid[plot_ping], Sv_38_regrid.range, label='Echolab2 Regrid', color='orange', linewidth=1)
plt.gca().invert_yaxis()
fig2.suptitle("Ping " + str(plot_ping) + " comparison EV vs Echolab2")
plt.xlabel("Sv (dB)")
plt.ylabel("Range (m)")
plt.legend()


fig2 = plt.figure()
plt.plot(Sv_38[plot_ping], Sv_38.range, label='Echolab Original', color='blue', linewidth=1)
plt.plot(Sv_38_regrid[plot_ping], Sv_38_regrid.range, label='Echolab2 Regrid', color='orange', linewidth=1)
plt.gca().invert_yaxis()
fig2.suptitle("Ping " + str(plot_ping) + " comparison Echolab2 vs Echolab2")
plt.xlabel("Sv (dB)")
plt.ylabel("Range (m)")
plt.legend()


fig2 = plt.figure()
plt.plot(ev_Sv_38[plot_ping], ev_Sv_38.range, label='EV', color='blue', linewidth=1)
plt.plot(ev_Sv_38_regrid[plot_ping], ev_Sv_38_regrid.range, label='EV Regrid', color='orange', linewidth=1)
plt.gca().invert_yaxis()
fig2.suptitle("Ping " + str(plot_ping) + " comparison EV vs EV")
plt.xlabel("Sv (dB)")
plt.ylabel("Range (m)")
plt.legend()


fig2 = plt.figure()
plt.plot(ev_Sv_38[plot_ping], ev_Sv_38.range, label='EV', color='blue', linewidth=1)
plt.plot(Sv_38[plot_ping], Sv_38.range, label='Echolab', color='orange', linewidth=1)
plt.gca().invert_yaxis()
fig2.suptitle("Ping " + str(plot_ping) + " comparison EV original vs Echolab original")
plt.xlabel("Sv (dB)")
plt.ylabel("Range (m)")
plt.legend()

plt.show()

print()
