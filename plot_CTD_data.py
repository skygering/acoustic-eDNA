###############################################
# Takes .evl files and plots depth vs date
# Cleans any file with depths above 0 to remove
# those points for plotting
# Skylar Gering 06/04/21
###############################################

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def plotDepth(x, y, filepath):
    """
    Plots depth over time from .evl files
    Input: x - Dates, y - depths, filepath - where to save the file
    Output: Saves matplotlib of y vs x to filepath
    """
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1) 
    plt.plot(x, y) # plot y vs x
    # Only have hour and minutes on x axis - samples usually taken on same day
    plt.gcf().autofmt_xdate()
    myFmt = mdates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(myFmt)

    plt.title('Depth over Time from ' + file)
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
bfile = inpath+'/'+'CTDtoEVL.list'
with open(bfile, 'r') as infile:
    CTD_files = infile.readlines() # read in each CDT file name and the toffset

for j in range(0, len(CTD_files)):
    # Get the filename and replace .asc with .evl - you now have all .evl filenames
    (file, _) = CTD_files[j].split(',')
    file = file.strip().replace('.asc', '.evl')
    infn = inpath+'/'+file
    
    # Read in each .evl file
    with open(infn, 'r', encoding = "ISO-8859-1") as infile:
        EVL_data = infile.readlines()

    # Delete file headers
    del EVL_data[0]
    del EVL_data[1]

    data = {}
    for line in EVL_data[1:]:
        line.rstrip()
        # Pull out the date, time, and depth (one additional variable ignored) -  all of these are strings
        (d, t, depth, tmpVars) = line.split(maxsplit=3)
        dt = datetime.strptime(d+t, "%Y%m%d%H%M%S%f") # use date and time strings to make datetime object
        data[dt] = -1*float(depth) # negate depth and make into float

    dt, depth = zip(*data.items()) # make date/time list and depth list from tuples 
    outfn = infn.replace('.evl', '.png')
    plotDepth(dt, depth, outfn)
    
    # If we have a max depth greater than 0, we have some bad data - this would be taking measurements above water
    # Remove the 'bad' points and plot without them
    if max(depth) > 0:
        # remove data that a depth greater than 0
        clean_dt,clean_depth= zip(*[(dt, depth) for dt, depth in data.items() if depth < 0])
        clean_outfn = inpath+'/CLEAN_'+file.replace('.evl', '.png')
        plotDepth(clean_dt, clean_depth, clean_outfn)

    


    
    

    
