using Dates
using Revise
using ElectroPhysiology, PhysiologyAnalysis

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose

fn = raw"G:\Data\Two Photon\2025-02-14-GRAB-DA\GRAB-DA2m-R1-R_da_puff_100um004.tif"
data2P = readImage(fn);
#deinterleave!(data2P) #This seperates the movies into two seperate movies
img_frames = get_all_frames(data2P)
img_stack = data2P.data_array[:,:,1]
#Calculate the stack dFoF
#moving average stack

import PhysiologyAnalysis.moving_average
dFoF = baseline_stack(img_stack)

dFoF |> size

using OffsetArrays
OffsetArrays.centered(ones(10,1))

#%%

time = data2P.t
og_y = mean(img_frames[:,:,:,1], dims = (1,2))[1,1,:]
dFoF = baseline_trace(og_y)

#%%
using Pkg; Pkg.activate("test")

using GLMakie, PhysiologyPlotting
fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, time, og_y, color = :blue)
ax2 = Axis(fig[2, 1])
lines!(ax2, time, dFoF, color = :magenta)
fig

#%% Now work on the ROI sections
#These parts are now from DeltaFF

#Go through all pixels and calculate the dFoF

