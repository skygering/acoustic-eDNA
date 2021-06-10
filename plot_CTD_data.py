###############################################
# Takes .evl files and plots depth vs date
# Cleans any file with depths above 0 to remove
# those points for plotting
# Skylar Gering 06/04/21
###############################################

from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

def evl_list(ctd_list):
    """
    From a CTDtoELV.list makes a list of .evl files that should be created 
    if SeaBird-wTime-toEchoviewLine.py was run
    Input: file path to CTDtoELV.list
    Output: list of corresponding .evl filenames
    """
    with open(ctd_list, 'r') as infile:
        ctd_files = infile.readlines() # read in each CDT file name and the toffset
    
    evl_files = []
    for ctd in ctd_files:
        # Get the filename and replace .asc with .evl - you now have all .evl filenames
        (file, _) = ctd.split(',')
        file = file.strip().replace('.asc', '.evl')
        infn = inpath+'/'+file
        evl_files.append(infn)
    return evl_files

def plotDepth(dt, depth, filepath):
    """
    Plots depth over time from .evl files
    Input: datetimes, depths, filepath - where to save the file
    Output: Saves matplotlib of y vs x to filepath
    """
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1) 
    plt.plot(dt, depth) # plot y vs x
    # Only have hour and minutes on x axis - samples usually taken on same day
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(myFmt)

    plt.title('Depth over Time from ' + os.path.basename(filepath))
    plt.xlabel('Date and Time (h:m)')
    plt.ylabel('Depth from surface (m)')
     
    # Change the number of x-axis ticks depending on how long the data timeframe is (less ticks for more time)
    min_int = 1
    if min(depth) < -350:
        min_int = 3
    elif min(depth) < -200:
        min_int = 2
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=min_int)) # to get a tick every min_int
    plt.savefig(filepath) # save figure
    plt.close()

### START HERE ###
inpath = '/Volumes/GeringSSD/GU201905_CTD' # read in data - CHANGE TO MATCH LOCAL FILE STRUCTURE
# Use CTDtoEVL list to plot - assumes that you just ran SeaBird-wTime-toEchoviewLine.py
ctd_list_file = inpath+'/'+'CTDtoEVL.list'

def main():
    evl_files = evl_list(ctd_list_file)
    for file in evl_files:
        # Read in each .evl file
        with open(file, 'r', encoding = "ISO-8859-1") as infile:
            evl_data = infile.readlines()

        depth_data = {}
        for line in evl_data[2:]: # exclude the first two lines which are headers
            line.rstrip()
            # Pull out the date, time, and depth (one additional variable ignored) -  all of these are strings
            (d, t, depth, tmpVars) = line.split(maxsplit=3)
            dt = datetime.strptime(d+t, "%Y%m%d%H%M%S%f") # use date and time strings to make datetime object
            depth_data[dt] = -1*float(depth) # negate depth and make into float

        dt, depth = zip(*depth_data.items()) # make date/time list and depth list from tuples 
        outfn = file.replace('.evl', ".png")
        plotDepth(dt, depth, outfn)

        # If we have a max depth greater than 0, we have some bad data - this would be taking measurements above water
        # Remove the 'bad' points and plot without them
        if max(depth) > 0:
            # remove data that a depth greater than 0
            clean_dt,clean_depth= zip(*[(dt, depth) for dt, depth in depth_data.items() if depth < 0])
            clean_outfn = file.replace('.evl', "") + "_CLEAN.png"
            plotDepth(clean_dt, clean_depth, clean_outfn)

if __name__ == "__main__":
    main()


    
    

    
