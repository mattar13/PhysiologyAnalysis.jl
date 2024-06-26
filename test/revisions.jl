using Revise
using Dates
using ElectroPhysiology 
using PhysiologyAnalysis
using Pkg; Pkg.activate("test")
#When we get to the dataframe creation we can fill in details here
using XLSX, DataFrames, Query
using FileIO, Images, ImageView, ImageMagick
import .PhysiologyAnalysis.traverse_root

#%% Open all dates for the imaging files
img_dir = raw"D:\Data\Calcium Imaging"
patch_dir = raw"D:\Data\Patching"
filename = raw"D:\Data\Analysis\data_analysis.xlsx"

#Don't often need to do this (eventually we need to add an update button)
#dataset = create2PDataSheet(img_dir, patch_dir; verbose = true)

PhysiologyAnalysis.__init__()
dataset = open2PDataSheet(filename)
all_datasheet = dataset["All Files"]
img_datasheet = dataset["TIF Files"]
patch_datasheet = dataset["ABF Files"] 
dataset["Paired Files"] = pair_experiments(patch_datasheet, img_datasheet)
save2PDataSheet(filename, dataset)


#%% Can we do an analysis for the IV curves
