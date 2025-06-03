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

experiment = readImage(img_da_10mM_fn)
stim_trace = readABF(stim_da_10mM_fn, flatten_episodic = true,stimulus_name = "IN 2", stimulus_threshold = 0.5)
stimulusProtocol = getStimulusProtocol(stim_trace)
addStimulus!(experiment, stimulusProtocol)
stim_ends = getStimulusEndTime(experiment)
dt = experiment.dt

#%%
fig = Figure(size = (1000, 500))

ax1a = Axis(fig[1,1], title = "Baseline with varying λ", xlabel = "Time (s)", ylabel = "F0")
ax2a = Axis(fig[2,1], xlabel = "Time (s)", ylabel = "dF/F")

#Calculate the mean trace
mean_trace = mean(experiment, dims = 1)[1,:,1]
F0 = mean_trace./mean(mean_trace)

# Plot original trace
lines!(ax1a, experiment.t, F0, color = :black, label = "Original")

# Plot baselines with different lambdas
baseline, dfof = baseline_trace(mean_trace, lam = 1e4, assym = 0.005, niter = 20)
lines!(ax1a, experiment.t, baseline, color = :red, label = "Drift")
lines!(ax2a, experiment.t, dfof, color = :red, label = "dF/F")

# Add legends to all axes
Legend(fig[1,2], ax1a, "Baselining", framevisible = false)
Legend(fig[2,2], ax2a, "dF/F Parameters", framevisible = false)

fig