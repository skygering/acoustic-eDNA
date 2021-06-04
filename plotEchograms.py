import numpy as np
from matplotlib.pyplot import figure, show, subplots_adjust, get_cmap
from echolab2.processing import processed_data, line
from echolab2.plotting.matplotlib import echogram
from echolab2.instruments import EK60

rawfiles = ['Volumes/GeringSSD/GU1905/EK60/GU19_05-D20191015-T112548.raw']

fig = figure()
# Set some properties for the sub plot layout.
subplots_adjust(left=0.11, bottom=0.1, right=0.98, top=.93, wspace=None,
                hspace=0.9)
