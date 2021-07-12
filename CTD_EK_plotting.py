from echolab2.processing import line
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from echolab2.plotting.matplotlib import echogram
import numpy as np

def plot_evl(ax, evl_infile, evl_path=""):
    '''
    plot_evl: plots the CTD track of one .evl file (depth vs. time) - does NOT plot acoustic data
    Inputs: ax (matplotlib pyplot axes object) - axes object from a plt.subplots() call
            evl_infile (string) - .evl infile name - if it don't have a path within the filename
                then the evl_path variable is needed
            png_outfile (string) - optional name of .png file to save with .png extension - if not provided picture isn't saved
            png_path (string) - optional variable required if png_outfile does not have a path
    Outputs:
            if show is True: plot will be displayed
    '''
    evl_line = line.read_evl(os.path.normpath(evl_path + "/" + evl_infile))
    dt = evl_line.ping_time
    depth = evl_line.data

    ax.plot(dt, depth)

    # formating date on x-axis
    plt.xticks(rotation=45)
    xfmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    # Change the number of x-axis ticks depending on how long the data timeframe is (less ticks for more time)
    min_int = 1
    if max(depth) > 350:
        min_int = 3
    elif max(depth) > 200:
        min_int = 2
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=min_int)) # to get a tick every min_int

    ax.set_title('Depth over Time from ' + os.path.basename(evl_infile))
    ax.set_xlabel("Time (H:M)")
    ax.set_ylabel("Depth (m)")
    ax.invert_yaxis()

def plot_echo(ax, ek, fq, fq_thresholds = [-90, -20], transducer_offset = 0.0):
    '''
    plot_echo: plots an echogram with a given EK80 object from pyEcholab, option of adding a CTD trace overlay, zooming
               in on the trace, and saving/showing the file
    Input: ek (EK80 object from pyEcholab) - ek object - must have already read in the raw file to the ek object
           ax (matplotlib.pyplot axes object) - axes object created from matplotlib.pyplot.subplots() call
           fq (integer) - frequency of the raw data to plot - 18000, 38000, 70000, 120000, 200000 are frequent options
           fq_thresholds (two element integer list) - [lower dB threshold, upper dB threshold] defines the range of decible
                                                    that are not considered noise
           transducer_offset (double) - optional transducer offset in meters
           show (boolean) - optional variable that determines if picture is shown, default is True
    Outputs: Returns the matplotlib pyplot axes object post plotting
             if show is True: plot will be displayed
             if png_outfile is provided the .png file will be saved to provided path
    '''
    ax.set_title(str(fq) + "Hz")
    raw_data_list = ek.get_channel_data(frequencies=fq)
    raw_data = raw_data_list[fq][0]
    for file in raw_data_list[fq][1:]: 
        # if there is more than one .raw file, this appends all data of specified frequency
        raw_data.append(file)

    cal_obj = raw_data.get_calibration()
    cal_obj.transducer_offset_z = transducer_offset
    Sv = raw_data.get_Sv(calibration=cal_obj, return_depth=True)
    echo_plot = echogram.Echogram(ax, Sv, threshold=[fq_thresholds[0],fq_thresholds[1]])
    return echo_plot

def plot_evl_trace(ax, echo_plot, trace_infn, trace_path = "", zoom = True):
    '''
    plot_evl_trace: add a evl depth trace to an echogram plot
    Inputs:
           echoplot - echogram object from pyEcholab
           ax (matplotlib.pyplot axes object) - axes object created from matplotlib.pyplot.subplots() call
           trace_infn (string path and filename) -  file name of evl CTD trace to overlay on echogram - this is optional
           trace_path (string) - optional variable required if trace_infn does not have a path
           zoom (boolean) - optional - if True will zoom in on echogram surrounding evl CTD trace
           show (boolean) - optional variable that determines if picture is shown, default is True
    Output: echo_plot - returns the updated echogram object
    '''
    trace_infn = os.path.normpath(trace_path + "/" + trace_infn)
    depth_line = line.read_evl(trace_infn)
    echo_plot.plot_line(depth_line, linewidth=2.5, color = "black")
    if zoom:
        echo_bottom = max(depth_line.data) * 1.25
        x_lims = (min(depth_line.ping_time).astype('float'), max(depth_line.ping_time).astype('float'))
        ax.set_ylim(echo_bottom, 0)
        ax.set_xlim(x_lims[0], x_lims[1])
    return echo_plot


def plot_segments(segments, title = "CTD Track: Depth vs Time"):
    '''
    plot_segments: plotting function designed for verification of segment seperation success after running 
                   create_segments_dic() or mark_usable_depth() 
    Inputs: segments (nested dictionary) - segment dictionary created by create_segments_dic function
            title (string) - optional string title for graph
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

    # plot aesthetics
    fig.autofmt_xdate()
    xfmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_title(title)
    ax.set_xlabel("Time (H:M)")
    ax.set_ylabel("Depth (m)")
    ax.invert_yaxis()
    plt.show()



