using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
using Pkg; Pkg.activate("test")
using GLMakie, PhysiologyPlotting
#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
println("Loading the quinpirole baseline data...")
img_da_10mM_fn = raw"H:\Data\Two Photon\2025-05-19-nirCAT-STR-DPUFF\nirCAT_s2_10mM_mq010.tif"
stim_da_10mM_fn = raw"H:\Data\Patching\2025-05-19-nirCAT-str\25519016.abf"

#Marla wants me to do a different baseline analysis. Lets work on some new ways to do this. Fitting an exponential to the baseline.

#We want to do an exponential and then subtract the baseline
#Open some data

experiment = readImage(img_da_10mM_fn)
stim_trace = readABF(stim_da_10mM_fn, flatten_episodic = true,stimulus_name = "IN 2", stimulus_threshold = 0.5)
stimulusProtocol = getStimulusProtocol(stim_trace)
addStimulus!(experiment, stimulusProtocol)
stim_ends = getStimulusEndTime(experiment)
dt = experiment.dt
#%%

fig = Figure(size = (1000, 500))

ax1a = Axis(fig[1,1], title = "With linear bridge", xlabel = "Time (s)", ylabel = "F0")
ax2a = Axis(fig[2,1], xlabel = "Time (s)", ylabel = "dF/F")
ax1b = Axis(fig[1,2], title = "With assymetric least squares", xlabel = "Time (s)", ylabel = "F0")
ax2b= Axis(fig[2,2], xlabel = "Time (s)", ylabel = "dF/F")

#Calculate the mean trace
mean_trace = mean(experiment, dims = 1)[1,:,1]
baseline_divisor = mean(mean_trace)
F0 = mean_trace ./ baseline_divisor
lines!(ax1a, experiment.t, F0, color = :black, label = "Original")
lines!(ax1b, experiment.t, F0, color = :black, label = "Original")

background = iterative_linear_bridge(experiment.t, F0, stim_ends)
lines!(ax1a, experiment.t, background, color = :red, label = "Linear Bridge")
lines!(ax2a, experiment.t, (F0 - background), color = :blue, label = "baseline subtracted")

background_als = PhysiologyAnalysis.baseline_als(F0, lam = 1e5, assym = 0.05, niter = 20)
baseline_trace_model = baseline_trace(mean_trace)
lines!(ax1b, experiment.t, background_als, color = :red, label = "ALS Baseline")
lines!(ax2b, experiment.t, baseline_trace_model, color = :blue, label = "baseline subtracted")


# Add legends to all axes
axislegend(ax1a, position = :rt)
axislegend(ax2a, position = :rt)
axislegend(ax1b, position = :rt)
axislegend(ax2b, position = :rt)

fig