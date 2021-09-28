# acoustic-eDNA
#### Summer project using acoustic and eDNA data to examine new fisheries monitoring technique as a NOAA Hollings Scholar

This project involves comparing the acoustic data and eDNA data collected on the 2019 Ecosystem Monitoring Survey (EcoMon). eDNA is a relatively new fisheries monitoring technique and while it has proven somewhat successful inland and in smaller water area, it is relativly untested in the open ocean and across the continental shelf and slope. We hope that by comparing the eDNA findings to the acoustic results we will be able to study to potential of eDNA monitoring and potentially provide a new form of acoustic data vertification.

This code base has two main python files: CTD_EK_processing.py and CTD_EK_plotting.py. You will also require the pyEcholab Package (https://github.com/CI-CMG/pyEcholab). If the RKT-80 branch is still live, use that one. If not, it has been merged into the main branch. 


**There are also three scripts that will walk you through the following process. These scripts, in order of use are `raw_overlay_ctd.py`, `segment_subsets_mfi.py`, and `mfi_masks_ABC.py`.** 

There are also other scripts in the test folder that are worth looking at for more guidance on how to apply the methods.

When starting, you need the .asc files output by the CTD trace and the .raw files that are output by the echosounder. This code is all based off of the EK80 echosounder software. Code will need to be updated if a different echosounder is used.

### Processing - CTD_EK_processing.py
First, you want to turn the .asc files into .evl files (files that could be read into/exported from Echoview). These just have the time and the depth of the CTD trace for each cast. For this, you also need a file called CTDtoEVL.list, which has the name of each .asc file followed by the time offset of the CTD readings from the echosounder readings. Here is an example from the first few lines of this file:

![ctd_to_evl_pic](https://user-images.githubusercontent.com/60117338/124821358-741f7480-df23-11eb-8446-cb249249544a.png)

To create the .evl files, run the `asc_to_evl` function in CTD_EK_processing.py. This function takes 4 parameters: the name of the file "CTDtoEVL.list", the path to both this list and the .asc files (they should be in the same place), the location to save the .evl files, and the Echoview version (there is already a default for this so unless you want to upload it to a specific Echoview version, don't worry about it). 

This will return files which have the following file structure:

![evl_example_pic](https://user-images.githubusercontent.com/60117338/124824773-ad59e380-df27-11eb-9d5a-261dd1ef5fee.png)

Where the first line is the Echoview version, the second is the number of points, and after that each line is the date of a point, the time, the depth, and the status (3 is good, 2 is bad).

If you want a list of the .asc files, you can run the  `asc_from_list function`, which also takes in the path to "CTDtoEVL.list". This will return a list of the .asc filenames without paths. This can be passed into the `cast_new_extension` function along with the current extension, desired extension, and an optional descriptor to add to the filename. To get a list of the .evl filenames, you can run:


`asc_names = asc_from_list("path/CTDtoEVL.list")`

`evl_names = cast_new_extension(asc_names, ".asc", ".evl")`

Now that we have the .evl files, we want to find the .raw files (the acoustic data) that overlaps with the CTD trace in the .evl files. For this, we will run the `match_raw_evl` function, which needs a list of all the .evl files and a list of all of the raw files (see script for a way to get all .raw files in a directory). This depends on the .raw file names being of the form: cruise-Dymd-Thms.raw. This function returns a dictionary with evl filenames as keys and a list of overlapping .raw files as values. If an outfile path is provided, this dictionary will be saved as a list called evl_raw_matches.list. This list can be read back into dictionary format with the function `evl_raw_dic_from_file`.  The output file should look something like the following:

![evl_raw_pic](https://user-images.githubusercontent.com/60117338/124824644-84395300-df27-11eb-83ea-bd7b6ac5720a.png)

We now want to split the CTD profile up into segments and identify at which segments eDNA samples were taken. Samples are only taken when the CTD trace plateau's at a single water depth for a prolonged time period. eDNA samples were taken at some of, but not all of, these plateaus. To analyze one cast, a list of depths where eDNA samples are taken is needed, as well as the .evl file for that cast. After analysis, a nested dictionary for the cast will be returned with points within the trace tagged as being in a water "bottle" sample and futher witin a "useable" eDNA sample. This dictionary can be saved as a nested .json file. While this can be done by hand by running first `create_segments_dic` and then `mark_usable_depth`, it is recommended that you use `interactive_segment_maker` as this will allow you to dynamically update the classification if it doesn't look correct. If you do use `interactive_segment_maker`, make sure to note the possible needed changes tothe `atol_zero` equation (explained in code comments). 

Before running any of this, you should make a file with each line being cast number followed by the depths at which an eDNA sample was taken. Use `read_casts` to read this file into a dictionary. Here are a few lines of an example file (if you want the method to work you need the header):

![edna_depth_pic](https://user-images.githubusercontent.com/60117338/124830175-5a375f00-df2e-11eb-8c60-1ea2d7b5ef2a.png)

Post running `interactive_segment_maker` you will have a .json file and/or a dictionary for each segment noting which points are in which segments and which segments are where water bottle samples and eDNA samples were taken. 

This segment file can be used to subset Sv data into smaller subsets that surround where the eDNA data was taken. The script `segment_subsets_mfi.py` creates 10 minute by 4 meter subsets. This can be done using the function `subset_segments_Sv`, which needs a segment dictionary as described above, as well as the .raw files that have the Sv data. This script allows you to easily create a box aound eDNA segments. It is recommended you use `interactive_subset_maker` as this allows you to dynamically adjust the bounds of the subset if it goes into the surface or the ocean floor. It also allows you to exclude some frequencies of data after seeing the noise. This outputs a dictionary, but it can be turned into a .json file using `subset_to_json`, which makes it easy to export for later use. 

Once you have the subsets, the next step is to calculate the MFI for each cast at each depth. This combines several frequencies of Sv data into a data array of biological classification. For more information see:
>Trenkel, Verena M., and Laurent Berger. "A fisheries acoustic multi-frequency indicator to inform on large scale spatial patterns of aquatic pelagic ecosystems."  

Once you have the MFI calculations, they can be used to create a mask on Sv data for the type of biological infomation you are interested in (0-0.4 swimbladder fish, 0.4-0.6 small resonant bubbles, 0.6-0.8 zooplankton, and 0.8-1 non-swimbladder fish). For this process, use `mask_mfi`. Finally, this masked Sv data can be used to calculate the area backscattering coefficent (ABC) using `calc_ABC`. 

We are planning on calculating the ABC for all usable casts and depths, seperating these values into quintiles, and then comaparing these to the 5 eDNA level rankings. 

### Plotting - CTD_EK_plotting.py

We are able to plot the CTD traces using `plot_evl`. For this, we just need a matplotlib pyplot axis to plot on and a .evl file.

We can also plot an echogram using `plot_echo`. This relies heavily on pyEcholab's `echogram` function and simply adds additional aesthetic features. 

We can also overlay a CTD trace over an echogram by passing the echogram returned by `plot_echo` into `plot_evl_trace` with an axis and a .evl trace.

We can plot a segment dictionary created by `create_segments_dic` and then `mark_usable_depth` in the same pattern used by the `interactive_segment_maker`. For this use `plot_segments`.

Finally, we can plot MFI seperated into the 4 biological catagories using `plot_MFI`.

### R Code

There is also some R code available to visualize the data. After running the python code to create the subset, MFI, or ABC data, these outputs can be read into R, and the availible R code ran to make the visualizations and do simple statistical analysis. These scripts were written to help anyone who strongly prefers working in R for analysis and visualization. 
