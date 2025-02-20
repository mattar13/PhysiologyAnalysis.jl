using Dates
using Revise
using ElectroPhysiology, PhysiologyAnalysis
using Statistics

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose

fn = raw"G:\Data\Two Photon\2025-02-14-GRAB-DA\GRAB-DA2m-R1-R_da_puff_100um004.tif"
ic_fn = raw"G:\Data\Patching\2025-02-14-da_puffs\25214001.abf"
data2P = readImage(fn);
addStimulus!(data2P, ic_fn, "IN 2")
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
idx_start = round(Int64, getStimulusEndTime(data2P)[idx]/data2P.dt)
idx_end = round(Int64, getStimulusStartTime(data2P)[idx+1]/data2P.dt)
segment_t = time2P[idx_start:idx_end]
segment_trace = df_trace[idx_start:idx_end]

fit = fit_parametric(segment_trace, segment_t)
y_fit = map(TIME-> single_stim_model(TIME, fit.param), segment_t)
fig, ax = lines(segment_t, segment_trace, color = :cyan)
lines!(ax, segment_t, y_fit, color = :black)

fit.resid