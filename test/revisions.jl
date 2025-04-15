using Dates
using Revise
using ElectroPhysiology, PhysiologyAnalysis
using Statistics

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose

fn = raw"E:\Data\Two Photon\2025-02-14-GRAB-DA\GRAB-DA2m-R1-R_da_puff_100um004.tif"
stim_fn = raw"E:\Data\Patching\2025-02-14-da_puffs\25214001.abf"

#fn = raw"F:\Data\Two Photon\2025-02-26-GRAB-DA\grab-da-R1-R-nirCAT-70um_kcl_AP-NOMF-GABA009.tif"
#stim_fn = raw"F:\Data\Patching\2025-02-26-da_puff\25226011.abf"

data2P = readImage(fn);
addStimulus!(data2P, stim_fn, "IN 2", flatten_episodic = true)
time2P = data2P.t

#deinterleave!(data2P) #This seperates the movies into two seperate movies
img_frames = get_all_frames(data2P)
img_stack = data2P.data_array[:,:,1]

dFoF = baseline_stack(data2P, window = 40)
df_trace = project(dFoF, dims = (1,2))[1,1,:,1]

#%%
using Pkg; Pkg.activate("test")
using GLMakie, PhysiologyPlotting
fig = Figure()

ax1 = Axis(fig[1, 1])
lines!(ax1, time2P, dFoF[200,:], color = :cyan)
lines!(ax1, time2P, df_trace, color = :black)
vlines!(ax1, getStimulusEndTime(data2P), color = :red)

fig

#%% Now work on the ROI fitting. 
#We want to pull out a single section from the stimulus

idx = 3
getStimulusEndTime(data2P)
idx_start = round(Int64, getStimulusEndTime(data2P)[idx]/data2P.dt)-50
idx_end = round(Int64, getStimulusStartTime(data2P)[idx+1]/data2P.dt)
segment_t = time2P[idx_start:idx_end]
segment_trace = dFoF[6000, idx_start:idx_end]

fit = fit_parametric(segment_trace, segment_t)
y_fit = map(TIME-> single_stim_model(TIME, fit.param), segment_t)
fig, ax = lines(segment_t, segment_trace, color = :cyan)
lines!(ax, segment_t, y_fit, color = :black)
fig

fit.param
fit.resid