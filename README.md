# acoustic-eDNA
#### Summer project using acoustic and eDNA data to examine new fisheries monitoring technique as a NOAA Hollings Scholar

This project involves comparing the acoustic data and eDNA data collected on the 2019 Gordon Gunther Cruise. eDNA is a relatively new fisheries monitoring technique and while it has proven somewhat successful inland and in smaller water area, it is relativly untested in the open water column. We hope that by comparing the eDNA findings to the acoustic results we will be able to study to potential effectivness of eDNA monitoring and potentiall improve practices.

This code base has two main python files: CTD_EK_processing.py and CTD_EK_plotting.py. You will also require the pyEcholab Package (https://github.com/CI-CMG/pyEcholab). If the RKT-80 branch is still live, use that one. If not, it has been merged into the main branch. There is also an example script (inert name here) that will walk you through the following process. 

When starting, you need the .asc files output by the CTD trace and the .raw files that are output by the echosounder. This code is all based off of the EK80 echosounder. Code may need to be updated if a different echosounder is used.

### Processing

First, you want to turn the .asc files into .evl files (files that could be read into/exported from Echoview). These just have the time and the depth of the CTD trace for each cast. For this, you also need a file called CTDtoEVL.list, which has the name of each .asc file followed by the time offset of the CTD readings from the echosounder readings. Here is an example from the first few lines of this file:

![ctd_to_evl_pic](https://user-images.githubusercontent.com/60117338/124821358-741f7480-df23-11eb-8446-cb249249544a.png)

To create the .evl files, run the `asc_to_evl` function in CTD_EK_processing.py. This function takes 4 parameters: the name of the file "CTDtoEVL.list", the path to both this list and the .asc files (they should be in the same place), the location to save the .evl files, and the Echoview version (there is already a default for this so unless you want to upload it to a specific Echoview version, don't worry about it). 

This will return files which have the following file structure:

![evl_example_pic](https://user-images.githubusercontent.com/60117338/124824773-ad59e380-df27-11eb-9d5a-261dd1ef5fee.png)

Where the first line is the Echoview version, the second is the number of points, and after that each line is the date of a point, the time, the depth, and the status (3 is good, 2 is bad).

If you want a list of the .asc files, you can run the  `asc_from_list function`, which also takes in the path to "CTDtoEVL.list". This will return a list of the .asc filenames without paths. This can be passed into the `cast_new_extension` function along with the current extension, desired extension, and an optional descriptor to add to the filename. To get a list of the .evl filenames, you can run:


`asc_names = asc_from_list("path/CTDtoEVL.list")`

`evl_names = cast_new_extension(asc_names, ".asc", ".evl")`

Now that we have the .evl files, we want to find the .raw files (the acoustic data) that overlaps with the CTD trace in the .evl files. For this, we will run the `match_raw_evl` function, which needs a list of all the .evl files and a list of all of the raw files (see script for a way to get all .raw files in a directory). This depends on the .raw file names being of the form: cruise-DYmd-Thms.raw. This function returns a dictionary with evl filenames as keys and a list of overlapping .raw files as values. If an outfile path is provided, this dictionary will be saved as a list called evl_raw_matches.list. This list can be read back into dictionary format with the function `evl_raw_dic_from_file`.  The output file should look something like the following:

![evl_raw_pic](https://user-images.githubusercontent.com/60117338/124824644-84395300-df27-11eb-83ea-bd7b6ac5720a.png)

We now want to split the CTD profile up into segments and identify at which segments eDNA samples were taken. Samples are only taken when the CTD trace plateau's at a single water depth for a prolonged time period. eDNA samples were taken at some of, but not all of, these plateaus. To analyze one case, a list of depths where eDNA samples are taken is needed, as well as the .evl file for that cast. After analysis, a nested dictionary for the cast will be returned with points within the trace tagged as being in a water "bottle" sample and futher witin a "useable" eDNA sample. This dictionary can be saved as a nested .json file. While this can be done by hand by running first `create_segments_dic` and then `mark_usable_depth`, it is recommended that you use `interactive_segment_maker` as this will allow you to dynamically update the classification if it doesn't look correct. If you do use `interactive_segment_maker`, make sure to note the possible needed changes tothe `atol_zero` equation (explained in code comments). 

Before running any of this, you should make a file with each line being cast number followed by the depths at which an eDNA sample was taken. Use `read_casts` to read this file into a dictionary. Here are a few lines of an example file (if you want the method to work you need the header):

![edna_depth_pic](https://user-images.githubusercontent.com/60117338/124830175-5a375f00-df2e-11eb-8c60-1ea2d7b5ef2a.png)

Post running `interactive_segment_maker` you will have a .json file and/or a dictionary for each segment noting which points are in which segments and which segments are where water bottle samples and eDNA samples were taken. 

Finally, this .json file or dictionary can be subset into 

### Plotting

Now that we have a dictionary of evl files to matching raw files, we can overlay the two to visualize them. 





