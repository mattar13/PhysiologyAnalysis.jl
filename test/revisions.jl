using Dates
using Revise
using ElectroPhysiology, PhysiologyAnalysis
using Statistics

using Pkg; Pkg.activate("test")
using GLMakie, PhysiologyPlotting

import ElectroPhysiology: Experiment, TWO_PHOTON
include("ROIVisualization.jl")

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
fn = raw"F:\Data\Two Photon\2025-05-02-GRAB-DA-nirCAT-STR\grab-nircat-str-kpuff_3x012.tif"
stim_fn = raw"F:\Data\Patching\2025-05-02-GRAB-DA-STR\25502017.abf"

data2P = readImage(fn);
addStimulus!(data2P, stim_fn, "IN 2", flatten_episodic = true)
deinterleave!(data2P) #This seperates the movies into two seperate movies
time2P = data2P.t

# Split the image into 8x8 pixel ROIs
pixel_splits_roi!(data2P, 8)

#%% Process all ROIs for channel 2 and stimulus 2
roi_analysis = process_rois(data2P; 
    channels=[1, 2],           # Only process channel 2
    stim_indices=[3],      # Only process the second stimulus
    delay_time=50.0,       # 50ms delay time for analysis
    sig_window=50.0,        # 50ms window to look for significant responses after stimulus
    window = 15,             # 15-point window for moving average
    n_stds = 5.0
)

# Store the analysis in the experiment's HeaderDict
# data2P.HeaderDict["ROI_Analysis"] = roi_analysis

# Get all significant ROIs and print summary
sig_rois = get_significant_rois(roi_analysis)
println("Found $(length(sig_rois)) significant ROIs")

# Get fit parameters and print summary statistics
fit_params = get_fit_parameters(roi_analysis)
println("Mean amplitude of significant ROIs: ", mean(first.(fit_params)))

fig = plot_roi_analysis(data2P, roi_analysis)
display(fig)

# Optionally save the figure
# save("roi_analysis.png", fig)
