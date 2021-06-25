from echolab2.processing import line
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def plot_evl(evl_infile, png_outfile = "", evl_path="", png_path="", show = True):
    '''
    plot_evl: plots the CTD track of one .evl file (depth vs. time)
    Inputs: evl_infile (string) - .evl infile name - if it don't have a path within the filename
                then the evl_path variable is needed
            png_outfile (string) - optional name of .png file to save with .png extension - if not provided picture isn't saved
            png_path (string) - optional variable required if png_outfile does not have a path
            show (boolean) - optional variable that determines if picture is shown, default is True
    Outputs:
            if show is True: plot will be displayed
            if png_outfile is provided the .png file will be saved to provided path
    '''
    evl_line = line.read_evl(os.path.normpath(evl_path + "/" + evl_infile))
    dt = evl_line.ping_time
    depth = evl_line.data

    fig, ax = plt.subplots()
    ax.plot(dt, depth)

    # formating date on x-axis
    fig.autofmt_xdate()
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

    if len(png_outfile) != 0:
        plt.savefig(os.path.normpath(png_path + "/" + png_outfile)) # save figure
    if show:
        plt.show()
    plt.close()