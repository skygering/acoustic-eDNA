import CTD_EK_processing as process
import CTD_EK_plotting as plotting
import find_plateau_CTD
import plot_CTD_data
import os

raw_path = "/Volumes/GeringSSD/GU1905_Acoustic/EK60"
ctd_path = '/Volumes/GeringSSD/GU201905_CTD'
ctd_list = '/CTDtoEVL.list'

# creates a list of the .asc CTD files names
with open(ctd_path + ctd_list, 'r') as infile:
    ctd_files = infile.readlines() # read in each CDT file name and the toffset
asc_files = []
for ctd in ctd_files:
    # Get the filename and replace .asc with .evl - you now have all .evl filenames
    (file, _) = ctd.split(',')
    asc_files.append(file)

# creates .evl files for each cast in CTDtoEVL.list and saves them in ctd_path
#process.asc_to_evl(ctd_list, ctd_path, ctd_path)

# create a list of the .evl file names, which were created in the above call
evl_list = process.cast_new_extension(asc_files, ".asc", ".evl")

# gets a list of the .raw files by giving a location and then finding all the .raw files in that directory
# this keeps the filepath in the filenames saved
raw_list = process.glob.glob(os.path.normpath(raw_path + "/*.raw"))

# finds overlapping .evl and .raw files and saves the list to a dictionary (also saves to file if given outfile_path)
#evl_raw_dic = process.match_raw_evl(evl_list, raw_list, evl_inpath="/Volumes/GeringSSD/GU201905_CTD")

png_list = process.cast_new_extension(evl_list, ".evl", ".png")
print(png_list)

plotting.plot_evl(evl_list[21], png_list[21], evl_path=ctd_path, png_path="/Volumes/GeringSSD")


eDNA_cast_dic = find_plateau_CTD.read_casts('/Volumes/GeringSSD/eDNA_cast.txt')
outfile_path = '/Volumes/GeringSSD/GU201905_CTD_JSON/'


evl_list = plot_CTD_data.evl_list(ctd_path+ctd_list)
for i in range(len(evl_list)):
    outfile = outfile_path + os.path.basename(evl_list[i]).replace(".evl", "") + "_segments.json"
    find_plateau_CTD.make_segments_json(eDNA_cast_dic[i+1], evl_list[i], outfile)