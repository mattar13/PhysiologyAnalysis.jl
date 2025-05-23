using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
#include("ROIVisualization.jl")

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
#We should look through the available files and see which ones fit
img_fn = raw"H:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b4_grabda-nircat-100uA_pulse012.tif"
stim_fn = raw"H:\Data\Patching\2025-05-15-GRAB-DA-STR\25515016.abf"

data = open2Pdata(img_fn, 
    ic_stim_filename = stim_fn, stimulus_name = "IN 3", 
    stimulus_threshold = 0.5, 
    split_channel = true, main_channel = :grn,
    spike_train = true,
    post_event_time = 120.0
)
xlims = data["xlims"]
ylims = data["ylims"]
experiment = data["experiment"]
dataIC = data["dataIC"]

# Split the image into 8x8 pixel ROIs
pixel_splits_roi!(experiment, 16)
# Process all ROIs for channel 2 and stimulus 2
roi_analysis = process_rois(experiment; 
    channels=[1, 2],           # Only process channel 2
    stim_indices=nothing,      # Only process the second stimulus
    delay_time=50.0,       # 50ms delay time for analysis
    sig_window=50.0,        # 50ms window to look for significant responses after stimulus
    window = 15,             # 15-point window for moving average
    n_stds = 5.0,
    lam = 1e4,
    assym = 0.075,
    niter = 100
)

roi_analysis |> typeof |> fieldnames
roi_analysis.analysis_parameters["delay_time"]
roi_analysis.analysis_parameters

# Extract significant ROIs and fit them
significant_rois = get_significant_rois(roi_analysis)
println("Found $(length(significant_rois)) significant ROIs")

process_significant_rois(roi_analysis)

dff = mean(hcat(get_dfof_traces(roi_analysis, significant_rois)...), dims =2 )[:, 1]
tseries = (1:length(dff)) .* experiment.dt
lines(tseries, dff)

roi_pixels = get_significant_roi_pixels(roi_analysis, experiment, channel_idx=2)

experiment