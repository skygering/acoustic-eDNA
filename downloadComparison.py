########################################################
# Compares data files gotten with wildcard terminal call (below call to get list of files followed by call to download files)
# > aws s3 ls s3://ncei-wcsd-archive/data/raw/Gordon_Gunter/GU1905/EK60/ --recursive --no-sign-request > lsfile
# > aws s3 cp s3://ncei-wcsd-archive/data/raw/Gordon_Gunter/ . --recursive --exclude "*" --include "GU1905*"  --no-sign-request
# to calls from Victoria's script gotten by mike
# Returns the difference between two lists of files
# Skylar Gering 06/04/21
########################################################

import numpy as np
from numpy.core.fromnumeric import size

# File that contains date, time, and size of download before the file names
filepath1 = "/Volumes/GeringSSD/mikeListOfFiles.txt"
with open(filepath1) as f:
    expected = [[token for token in line.split()] for line in f.readlines()] # read each line
    expectedArray = expected[7:-3] # first 7 and last 3 lines are headers/footers
    expectedFiles = set() # create a set of expected files
    for i in range(len(expectedArray)):
        expectedFiles.add(expectedArray[i][3]) # add each file name from mike's list

    with open("/Volumes/GeringSSD/listOfFiles.txt") as f1: # file with just list of files
        acutalFiles = set()
        acutalFiles.update([line.strip() for line in f1.readlines()]) # read in all lines and add to set
        
        diff = list(expectedFiles - acutalFiles) + list(acutalFiles - expectedFiles) # check difference between sets
        print(diff) # print difference - should be 0 or non .raw files
