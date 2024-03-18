using ElectroPhysiology, PhysiologyAnalysis
using Pkg; Pkg.activate("test")
using GLMakie
using FileIO, ImageView, Images
using Statistics
#%% Point to the file location ____________________________________________________________#
root = "F:/Data/Calcium Images/"
folder = "2024_03_11_VglutGC6/Cell3/"
file = "freeRec001.tif"
filepath = joinpath(root, folder, file);
filepath
# Open and prepare the data______________________________________________________________#
fps = 1.45
data = readImage(filepath; sampling_rate = fps);
px_x, px_y = data.HeaderDict["framesize"]
xlims = data.HeaderDict["xrng"]; #Extract the x domain
ylims = data.HeaderDict["yrng"]; #Extract the y domain
t = data.t; #Extract the t domain

mov = get_all_frames(data); #get all frames as a 3D array
fluo = mean(mov, dims = (1,2))[1,1,:]
zproj = mean(mov, dims = 3)[:,:,1] #Can we save this directly? 