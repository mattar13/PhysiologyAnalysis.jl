using Dates
using Revise
using ElectroPhysiology, PhysiologyAnalysis
using Statistics

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose

fn = raw"F:\Data\Two Photon\2025-02-14-GRAB-DA\GRAB-DA2m-R1-R_da_puff_100um004.tif"
data2P = readImage(fn);
ic_fn = raw"F:\Data\Patching\2025-02-14-da_puffs\25214001.abf"
dataIC = readABF(ic_fn, flatten_episodic = true) #Open the IC data

#Set the time offset
start2P = data2P.HeaderDict["FileStartDateTime"]-Second(3.0) #The computer clocks are off by 3 seconds
startIC = dataIC.HeaderDict["FileStartDateTime"]
t_offset = Millisecond(startIC - start2P).value/1000 #This is the delay between me pressing the buttons
dataIC.t .+= t_offset

#deinterleave!(data2P) #This seperates the movies into two seperate movies
img_frames = get_all_frames(data2P)
img_stack = data2P.data_array[:,:,1]

#%% Calculate the stack dFoF
import PhysiologyAnalysis.baseline_stack

dFoF = baseline_stack(data2P)
df_trace = project(dFoF, dims = (1,2))[1,1,:,1]

time2P = data2P.t
timeIC = dataIC.t
ic_trace = dataIC.data_array[1,:,3]

xcorr_xy = cor_xy(data2P.t, dF_trace, dataIC.t, dataIC.data_array[1,:,3])

#%%
using Pkg; Pkg.activate("test")
using GLMakie, PhysiologyPlotting
fig = Figure()

ax1 = Axis(fig[1, 1])
lines!(ax1, time2P, df_trace, color = :blue)

ax2 = Axis(fig[2, 1])
lines!(ax2, timeIC, ic_trace, color = :magenta)
linkxaxes!(ax1, ax2)

ax3 = Axis(fig[1, 2])
lines!(ax3, xcorr_xy, color = :blue)

fig
#%% Now work on the ROI sections
#These parts are now from DeltaFF

#Go through all pixels and calculate the dFoF

