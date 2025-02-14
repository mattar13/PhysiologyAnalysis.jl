using Dates
using Revise
using ElectroPhysiology, PhysiologyAnalysis

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose

fn = raw"G:\Data\Two Photon\2025-02-12-GRAB-DA\grab-da-r2_da-puff10mM014.tif"
data2P = readImage(fn);
#deinterleave!(data2P) #This seperates the movies into two seperate movies

img_frames = get_all_frames(data2P)
pixel_splits(size(img_frames)[1:2], 15)

time = data2P.t
og_y = mean(img_frames[:,:,:,1], dims = (1,2))[1,1,:]
baseline_y = baseline_als(og_y)
include("../src/Analysis/ImagingAnalysis/NanoImg.jl")
baseline_ma_y = baseline_ma(og_y-baseline_y)

#%%
using Pkg; Pkg.activate("test")

using GLMakie, PhysiologyPlotting
fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, time, og_y, color = :blue)
lines!(ax1, time, baseline_y, color = :red)
lines!(ax1, time, baseline_ma_y, color = :green)
ax2 = Axis(fig[2, 1])
lines!(ax1, time, og_y - baseline_y, color = :magenta)
fig