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

#%% Make plots with the range of baseline parameters
fig = Figure(size = (1800, 600))

# Lambda parameter plots
gl_lambda = fig[1, 1] = GridLayout()
ax1a = Axis(gl_lambda[1,1], title = "Baseline with varying λ", xlabel = "Time (s)", ylabel = "F0")
ax2a = Axis(gl_lambda[2,1], xlabel = "Time (s)", ylabel = "dF/F")

# Asymmetry parameter plots
gl_p = fig[1, 2] = GridLayout()
ax1b = Axis(gl_p[1,1], title = "Baseline with varying p", xlabel = "Time (s)", ylabel = "F0")
ax2b = Axis(gl_p[2,1], xlabel = "Time (s)", ylabel = "dF/F")

# Iterations parameter plots
gl_niter = fig[1, 3] = GridLayout()
ax1c = Axis(gl_niter[1,1], title = "Baseline with varying iterations", xlabel = "Time (s)", ylabel = "F0")
ax2c = Axis(gl_niter[2,1], xlabel = "Time (s)", ylabel = "dF/F")

#Calculate the mean trace
mean_trace = mean(experiment, dims = 1)[1,:,1]
F0 = mean_trace./mean(mean_trace)

# Create parameter ranges
lambdas = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
p_values = [0.001, 0.01, 0.05, 0.1, 0.25, 0.5]
niter_values = [5, 10, 20, 50, 100, 200]

# Create color gradients
colors_lambda = cgrad(:viridis, length(lambdas), categorical = true);
colors_p = cgrad(:plasma, length(p_values), categorical = true);
colors_niter = cgrad(:magma, length(niter_values), categorical = true);

# Plot original traces
lines!(ax1a, experiment.t, F0, color = :black, label = "Original")
lines!(ax1b, experiment.t, F0, color = :black, label = "Original")
lines!(ax1c, experiment.t, F0, color = :black, label = "Original")

# Plot with varying lambda
for (i, λ) in enumerate(lambdas)
    baseline, dfof = baseline_trace(mean_trace, lam = λ, assym = 0.05, niter = 20)
    lines!(ax1a, experiment.t, baseline, color = colors_lambda[i], label = "λ = $(Int(λ))")
    lines!(ax2a, experiment.t, dfof, color = colors_lambda[i], label = "λ = $(Int(λ))")
end

# Plot with varying p
for (i, p) in enumerate(p_values)
    baseline, dfof = baseline_trace(mean_trace, lam = 1e5, assym = p, niter = 20)
    lines!(ax1b, experiment.t, baseline, color = colors_p[i], label = "p = $p")
    lines!(ax2b, experiment.t, dfof, color = colors_p[i], label = "p = $p")
end

# Plot with varying niter
for (i, niter) in enumerate(niter_values)
    baseline, dfof = baseline_trace(mean_trace, lam = 1e5, assym = 0.05, niter = niter)
    lines!(ax1c, experiment.t, baseline, color = colors_niter[i], label = "niter = $niter")
    lines!(ax2c, experiment.t, dfof, color = colors_niter[i], label = "niter = $niter")
end

# Create separate legends
Legend(gl_lambda[1:2, 2], ax1a, "Lambda Parameters", framevisible = false)
Legend(gl_p[1:2, 2], ax1b, "Asymmetry Parameters", framevisible = false)
Legend(gl_niter[1:2, 2], ax1c, "Iteration Parameters", framevisible = false)

fig