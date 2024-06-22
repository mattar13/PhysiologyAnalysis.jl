using Revise
using ElectroPhysiology 
using PhysiologyAnalysis
using Pkg; Pkg.activate("test")
#When we get to the dataframe creation we can fill in details here
using XLSX, DataFrames, Query
using FileIO, Images, ImageView, ImageMagick
import .PhysiologyAnalysis.traverse_root

patch_dir = raw"D:\Data\Patching"
CaImg_dir = raw"D:\Data\Calcium Imaging"
patch_files = traverse_root(patch_dir)
twoP_files = traverse_root(CaImg_dir)
foreach(println, patch_files)
PhysiologyAnalysis.__init__()


date_test = readABF(patch_files[1])
img_test = readImage(twoP_files[1])