using Dates
using Revise
using ElectroPhysiology, PhysiologyAnalysis

include("../src/Analysis/ImagingAnalysis/NanoImg.jl")
#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose

fn = raw"G:\Data\Two Photon\2025_01_08_VGGC6_SWCNT\swcnt2_r3_1-50bo_kpuff_nomf_rac_10um011.tif"
data2P = readImage(fn);
deinterleave!(data2P) #This seperates the movies into two seperate movies

img_frames = get_all_frames(data2P)
pixel_splits(size(img_frames)[1:2], 15)

time = data2P.t
og_y = img_frames[100,100,:,2]
baseline_y = baseline_als(og_y)

#%%
using Pkg; Pkg.activate("test")

using GLMakie, PhysiologyPlotting
fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, time, og_y, color = :blue)
lines!(ax1, time, baseline_y, color = :red)