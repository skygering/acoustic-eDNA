from echolab2.processing import line
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from echolab2.plotting.matplotlib import echogram
import numpy as np
import CTD_EK_processing as process
import math
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker

def plot_evl(ax, evl_infile, evl_path="", title = ""):
    '''
    plot_evl: plots the CTD track of one .evl file (depth vs. time) - does NOT plot acoustic data
    Inputs: ax (matplotlib pyplot axes object) - axes object from a plt.subplots() call
            evl_infile (string) - .evl infile name - if it don't have a path within the filename
                then the evl_path variable is needed
            png_outfile (string) - optional name of .png file to save with .png extension - if not provided picture isn't saved
            png_path (string) - optional variable required if png_outfile does not have a path
            title (string) - title to be displayed at the top of plot - default is evl_infile: CTD Profile
    Outputs:
            After running function plt.show(), plt.savefig(), or plt.close() can all be run
    '''
    print("Plotting: " + evl_infile)
    evl_line = line.read_evl(os.path.normpath(evl_path + "/" + evl_infile))
    dt = evl_line.ping_time
    depth = evl_line.data

    ax.plot(dt, depth)
    # title
    if len(title) ==0:
        title = os.path.basename(evl_infile) + ": CTD Profile"
    ax.set_title(title)
    ax.set_xlabel("Time (H:M)")
    ax.set_ylabel("Depth (m)")
    # formating date on x-axis
    ax.tick_params(axis='x', labelrotation=15)
    xfmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    # Change the number of x-axis ticks depending on how long the data timeframe is (less ticks for more time)
    min_int = 2
    if max(depth) > 350:
        min_int = 3
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=min_int)) # to get a tick every min_int
    ax.invert_yaxis()

def plot_echo(ax, ek, fq, fq_thresholds = [-90, -20], transducer_offset = 0.0, title = ""):
    '''
    plot_echo: plots an echogram with a given EK80 object from pyEcholab, option of adding a CTD trace overlay, zooming
               in on the trace, and saving/showing the file
    Input: ek (EK80 object from pyEcholab) - ek object - must have already read in the raw file to the ek object
           ax (matplotlib.pyplot axes object) - axes object created from matplotlib.pyplot.subplots() call
           fq (integer) - frequency of the raw data to plot - 18000, 38000, 70000, 120000, 200000 are frequent options
           fq_thresholds (two element integer list) - [lower dB threshold, upper dB threshold] defines the range of decible
                                                    that are not considered noise
           transducer_offset (double) - optional transducer offset in meters
           title (string) - title to be displayed at the top of plot - default is frequency Hz
    Outputs: Returns echogram object - plot_evl-trace takes in this object
             After running function plt.show(), plt.savefig(), or plt.close() can all be run
    '''
    print("Plotting: echogram " + str(fq) + "Hz")
    if len(title) == 0:
        title = str(fq) + "Hz"
    ax.set_title(title)
    ax.tick_params(axis='x', labelrotation=15)
   
    Sv, _ = process.raw_to_Sv(ek, fq, transducer_offset)
    echo_plot = echogram.Echogram(ax, Sv, threshold=[fq_thresholds[0],fq_thresholds[1]])
    return echo_plot

def plot_evl_trace(ax, echo_plot, trace_infn, trace_path = "", zoom = True, time_offset = [2, 0]):
    '''
    plot_evl_trace: add a evl depth trace to an echogram plot
    Inputs:
           echo_plot - echogram object from pyEcholab
           ax (matplotlib.pyplot axes object) - axes object created from matplotlib.pyplot.subplots() call
           trace_infn (string path and filename) -  file name of evl CTD trace to overlay on echogram - this is optional
           trace_path (string) - optional variable required if trace_infn does not have a path
           zoom (boolean) - optional - if True will zoom in on echogram surrounding evl CTD trace
           show (boolean) - optional variable that determines if picture is shown, default is True
    Output: echo_plot - returns the updated echogram object
            After running function plt.show(), plt.savefig(), or plt.close() can all be run
    '''
    print("Adding CTD profile: " + trace_infn)
    trace_infn = os.path.normpath(trace_path + "/" + trace_infn)
    depth_line = line.read_evl(trace_infn)
    echo_plot.plot_line(depth_line, linewidth=2.5, color = "black")
    if zoom:
        echo_bottom = max(depth_line.data) * 1.35
        x_lims = ((min(depth_line.ping_time) - np.timedelta64(time_offset[0], 'm')).astype('float'), 
                 (max(depth_line.ping_time) + np.timedelta64(time_offset[1], 'm')).astype('float'))
        ax.set_ylim(echo_bottom, 0)
        ax.set_xlim(x_lims[0], x_lims[1])
    return echo_plot


def plot_segments(segments, title = "CTD Profile Segments", show = True):
    '''
    plot_segments: plotting function designed for verification of segment seperation success after running 
                   create_segments_dic() or mark_usable_depth() 
    Inputs: segments (nested dictionary) - segment dictionary created by create_segments_dic function
            title (string) - optional string title for graph
            show (boolean) - if True, show segments - needed for interactive segment finder
    Outputs: shows a graph with ascents and decents marked in red, plateaus in blue, and usable segments
             with dotted blue for easy vertification
    Note: This is used in CTD_EK_processing.py for segment creation vertification
    '''
    fig, ax = plt.subplots()
    for num in segments:
        style = 'r-'
        if segments[num]["bottle"]:
            style = 'b-'
            if segments[num]["usable"]:
                style = 'b:'
        x_seg, y_seg = zip(*segments[num]["points"])
        x_seg = [np.datetime64(date) for date in x_seg]
        ax.plot(x_seg, y_seg, style)
        ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=3))

    # plot aesthetics
    fig.autofmt_xdate()
    xfmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_title(title)
    ax.set_xlabel("Time (H:M)")
    ax.set_ylabel("Depth (m)")
    ax.invert_yaxis()
    ax.tick_params(axis='x', labelrotation=15)
    if show:
        plt.show()

def plot_MFI(ax, mfi, title = "", label_size = 14):
    '''
    plot_MFI: plot MFI processed data object created from calc_MFI()
    Inputs: ax (matplotlib.pyplot axes object) - axes object created from matplotlib.pyplot.subplots() call
            mfi (MFI processed data object) - see calc_MFI() and Rick Towler's processed_data.py
            title (string) - title of the plot
            label_size (integer) - size of label for color bar
    Outputs: image object created from imshow()
    Note: based off of Rick Towler's echogram.Echogram()
    '''
    def format_datetime(x, pos=None):
        '''
        format_datetime: attempts to convert floats into datetime64 objectss
        '''
        try:
            dt = x.astype('datetime64[ms]').astype('object')
            tick_label = dt.strftime("%H:%M:%S")
        except:
            tick_label = ''
        return tick_label

    # splits colors into 4 MFI catagories and value ranges
    colors = ["#006164", "#57c4ad", "#eda247", "#db4325"]
    bounds = [0, 0.4, 0.6, 0.8, 1]
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(bounds, 4)
    cmap.set_bad(color="gray")

    # rotates MFI data to be plotted depth vs time
    mfi_data = np.flipud(np.rot90(mfi.data, 1))

    # axis ticks
    yticks = mfi.depth
    xticks = mfi.ping_time.astype('float')
    # plot
    mfi_image = ax.imshow(mfi_data, cmap=cmap, norm = norm, aspect='auto', interpolation='none', 
                extent=[xticks[0], xticks[-1], yticks[-1], yticks[0]], origin='upper')

    # axis aesthetics
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_datetime))
    ax.tick_params(axis='x', labelrotation=15)
    y_label = 'Depth (m)'
    try:
        x = ax.get_xticks()[0]
        dt = x.astype('datetime64[ms]').astype('object')

        x_label = dt.strftime("%m-%d-%Y")
    except:
        x_label = ''
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(True, color='k')
    #title
    if len(title) == 0:
        title = "MFI Predictions"
    ax.set_title(title)
    # color bar aesthetics
    cbar = plt.colorbar(mfi_image)
    cbar.set_ticks(list())
    for index, label in enumerate(["Swimbladder Fish", "Small Resonant Bubbles", "Zooplankton", "Non-Swimbladder Fish"]):
        x = 1.5
        y = (2*index+1)/8
        cbar.ax.text(x,y,label, size = label_size)
    return mfi_image
