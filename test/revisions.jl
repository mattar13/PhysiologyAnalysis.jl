using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
#include("ROIVisualization.jl")

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
#We should look through the available files and see which ones fit
img_fn = raw"F:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b3_grabda-nircat-100uA_pulse_NOMF010.tif"
stim_fn = raw"F:\Data\Patching\2025-05-15-GRAB-DA-STR\25515012.abf"
println("Loading the nomf data")
data_nomf = open2Pdata(img_fn, 
    ic_stim_filename = stim_fn, stimulus_name = "IN 3", 
    spike_train = true,
    split_channel = true, main_channel = :grn,
    stimulus_threshold = 0.5, 
    post_event_time = 120.0
)
xlims = data_nomf["xlims"]
ylims = data_nomf["ylims"]
exp_nomf = data_nomf["experiment"]
stim_nomf = data_nomf["dataIC"]

pixel_splits_roi!(exp_nomf, 16)
roi_analysis = process_rois(exp_nomf; 
    channels=[1, 2],       # Only process channel 2
    stim_indices=nothing,  # Only process the second stimulus
    delay_time=50.0,       # 50ms delay time for analysis
    sig_window=50.0,       # 50ms window to look for sig resp after stimulus
    window = 15,           # 15-point window for moving average
    n_stds = 5.0, 
    lam = 1e4,  #These are baselineing parameters
    niter = 100
)

# Extract significant ROIs and fit them
significant_rois = get_significant_rois(roi_analysis; channel_idx=2)
println("Found $(length(significant_rois)) significant ROIs")

dff = mean(hcat(get_dfof_traces(roi_analysis, significant_rois; channel_idx=2)...), dims =2 )[:, 1]
tseries = (1:length(dff)) .* exp_nomf.dt