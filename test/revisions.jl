using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
using Pkg; Pkg.activate("test")
using GLMakie

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
println("Loading the quinpirole baseline data...")
img_fn = raw"F:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b6_grabda-nircat-300uA_pulse016.tif"
stim_fn = raw"F:\Data\Patching\2025-05-15-GRAB-DA-STR\25515025.abf"

data_quin_base = load_electric_data(img_fn, stim_fn, red_lam = 1e3, grn_lam = 1e3)

all_rois = data_quin_base["sig_traces"][2][3]
lines(all_rois[1])
vlines!([50.0], color = :red)
for i in 2:length(all_rois)
    lines!(all_rois[i])
end


#%% create a multidimensional array
n_rois = length(data_quin_base["sig_traces"][1][1])
n_stims = length(data_quin_base["sig_traces"][1])
n_timepoints = length(data_quin_base["sig_traces"][1][1][1])
n_channels = length(data_quin_base["sig_traces"])

sig_traces = zeros(n_stims, n_rois, n_timepoints, n_channels)

for STIM in 1:n_stims
    for CHANNEL in 1:n_channels
        for ROI in 1:n_rois
            println("STIM: $STIM, CHANNEL: $CHANNEL, ROI: $ROI")
            trace = hcat(data_quin_base["sig_traces"][CHANNEL][STIM][ROI]...)
            println(trace |> size)
            #find an NaNs
            nan_idx = findfirst(isnan.(trace))
            println("NaN at $nan_idx")
            sig_traces[STIM, ROI, :, CHANNEL] = trace
        end
    end
end


mean_trace = mean(sig_traces, dims = (1, 2))[:,1,:,:]
lines(mean_trace[1,:,2])

lines(data_quin_base["dff_red_trace"])
vlines!(data_quin_base["pks"], color = :red)

#%%
sig_traces = hcat(data_quin_base["sig_traces"][1][3]...)'
data_quin_base["sig_tseries"][1]

roi_analysis = data_quin_base["roi_analysis"]
sig_rois = get_significant_rois(roi_analysis, channel_idx = 1)

dfof_traces = get_dfof_traces(roi_analysis, sig_rois, stim_idx = 1, channel_idx = 1)

