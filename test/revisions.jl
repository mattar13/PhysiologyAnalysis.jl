using Dates
using Revise
using ElectroPhysiology, PhysiologyAnalysis
using Statistics

using Pkg; Pkg.activate("test")
using GLMakie, PhysiologyPlotting

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
fn = raw"H:\Data\Two Photon\2025-05-02-GRAB-DA-nirCAT-STR\grab-nircat-str-kpuff_3x012.tif"
stim_fn = raw"H:\Data\Patching\2025-05-02-GRAB-DA-STR\25502017.abf"

data2P = readImage(fn);
addStimulus!(data2P, stim_fn, "IN 2", flatten_episodic = true)
deinterleave!(data2P) #This seperates the movies into two seperate movies
time2P = data2P.t
pixel_splits_roi!(data2P, 8)

data_dict = roi_processing(data2P, 2, channel_index = 1)
data_dict[2]["dFoF"]
data_dict[2]["fit_param"]
data_dict[2]["significance"]

t_series = data_dict[2]["t_series"]

#%% Find a good way to plot the data
sigs = get_significance(data_dict)
dfof = get_dFoF(data_dict) 

findall(sigs .!= 0.0)
sig_traces = mean(dfof[findall(sigs .!= 0.0), :], dims = 1)[1,:]

#%% Makes some plots
fig = Figure()
ax = Axis(fig[1, 1], title = "dFoF", xlabel = "Time (s)", ylabel = "dFoF")
lines!(ax, t_series, sig_traces, color = :blue, label = "dFoF")

fig