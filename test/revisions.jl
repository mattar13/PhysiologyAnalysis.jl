using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
#include("ROIVisualization.jl")

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
println("Loading the quinpirole baseline data...")
img_fn = raw"H:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b6_grabda-nircat-300uA_pulse016.tif"
stim_fn = raw"H:\Data\Patching\2025-05-15-GRAB-DA-STR\25515025.abf"

data_quin_base = load_electric_data(img_fn, stim_fn)

#%% create a multidimensional array
n_rois = length(data_quin_base["sig_traces"][1][1])
n_stims = length(data_quin_base["sig_traces"][1])
n_timepoints = length(data_quin_base["sig_traces"][1][1][1])
n_channels = length(data_quin_base["sig_traces"])

sig_traces = fill(nothing, n_stims, n_timepoints, n_channels)

for STIM in 1:n_stims
    for CHANNEL in 1:n_channels
        println("STIM: $STIM, CHANNEL: $CHANNEL")
        trace = hcat(data_quin_base["sig_traces"][CHANNEL][STIM]...)
        println(trace |> size)
        # sig_traces[i,j,k,l] = data_quin_base["sig_traces"][l][i][j][k]
    end
end

sig_traces

#%%
sig_traces = hcat(data_quin_base["sig_traces"][1][3]...)'
data_quin_base["sig_tseries"][1]

roi_analysis = data_quin_base["roi_analysis"]
sig_rois = get_significant_rois(roi_analysis, channel_idx = 1)

dfof_traces = get_dfof_traces(roi_analysis, sig_rois, stim_idx = 1, channel_idx = 1)

