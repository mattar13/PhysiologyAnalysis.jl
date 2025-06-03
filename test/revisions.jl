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

data = open2Pdata(img_da_10mM_fn, 
    #Indicate the stimulus dataset
    stim_filename = stim_da_10mM_fn, #New option if we want to specify a 2P stimulus dataset
    stimulus_name = "IN 2", 
    stimulus_threshold = 0.5,
    
    split_channel = false, 
    main_channel = :red,
    pre_event_time = 20.0, 
    post_event_time = 60.0, 

    verbose = 3, # 1: basic progress, 2: timing info, 3: detailed info
)

#%% ╔═╡Now make example plots to show the results
fig = Figure(size = (1000, 500))

ax1a = Axis(fig[1,1], title = "Original", xlabel = "Time (s)", ylabel = "dF/F")
ax2a = Axis(fig[2,1], title = "dF/F", xlabel = "Time (s)", ylabel = "dF/F")

ax1b = Axis(fig[1,2], title = "Original", xlabel = "Time (s)", ylabel = "dF/F")
ax2b = Axis(fig[2,2], title = "dF/F", xlabel = "Time (s)", ylabel = "dF/F")

lines!(ax1a, data["time"], data["grn_trace"], color = :green, label = "Green")
lines!(ax2a, data["time"], data["dff_grn_trace"], color = :green, label = "Green")

lines!(ax1b, data["time"], data["red_trace"], color = :red, label = "Red")
lines!(ax2b, data["time"], data["dff_red_trace"], color = :red, label = "Red")

fig

#%% Now do the ROI analysis
data = PhysiologyAnalysis.load_puffing_data(img_da_10mM_fn, stim_da_10mM_fn)