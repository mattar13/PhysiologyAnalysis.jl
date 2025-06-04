using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
using Pkg; Pkg.activate("test")
using GLMakie, PhysiologyPlotting

#%% ╔═╡Found an issue. With really crazy spikes, the baseline correction is not working.
println("Loading the quinpirole baseline data...")
img_fn = raw"H:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b4_grabda-nircat-100uA_pulse012.tif"
stim_fn = raw"H:\Data\Patching\2025-05-15-GRAB-DA-STR\25515016.abf"

# img_fn = raw"H:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b6_grabda-nircat-300uA_pulse016.tif"
# stim_fn = raw"H:\Data\Patching\2025-05-15-GRAB-DA-STR\25515025.abf"

data = open2Pdata(img_fn, 
    #Indicate the stimulus dataset
    stim_filename = stim_fn, #New option if we want to specify a 2P stimulus dataset
    stimulus_name = "IN 3", 
    spike_train = true,
    
    split_channel = true, 
    main_channel = :red,
    pre_event_time = 20.0, 
    post_event_time = 60.0, 
    red_spike_reduction = :median,
    grn_spike_reduction = :median,
    red_window = 100,
    grn_window = 5,
    red_lam = 1e4,
    red_assym = 0.005,
    verbose = 3, # 1: basic progress, 2: timing info, 3: detailed info
)
# ╔═╡Now make example plots to show the results
fig = Figure(size = (1000, 500))

ax1a = Axis(fig[1,1], title = "Original", xlabel = "Time (s)", ylabel = "dF/F")
ax2a = Axis(fig[2,1], title = "dF/F", xlabel = "Time (s)", ylabel = "dF/F")

lines!(ax1a, data["time"], data["red_trace"]./mean(data["red_trace"]), color = :red, label = "Red")
lines!(ax1a, data["time"], data["red_drift"], color = :blue, label = "Drift")
vlines!(ax1a, data["tstamps"], color = :black, linestyle = :dash)

lines!(ax2a, data["time"], data["dff_red_trace"], color = :green, label = "dF/F")
vlines!(ax2a, data["tstamps"], color = :black, linestyle = :dash)

fig