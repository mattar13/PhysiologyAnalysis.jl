using Revise
using ElectroPhysiology 
using PhysiologyAnalysis
using Pkg; Pkg.activate("test")
#When we get to the dataframe creation we can fill in details here
using XLSX, DataFrames, Query

import .PhysiologyAnalysis.traverse_root

patch_dir = raw"D:\Data\Patching"
CaImg_dir = raw"D:\Data\Calcium Imaging"
patch_files = traverse_root(patch_dir)
twoP_files = traverse_root(CaImg_dir)

#I want to extract the date from .abf files
data_test = ElectroPhysiology.readABFInfo(patch_files[1])
#Now we pull out FileStartDateTime for each .abf file