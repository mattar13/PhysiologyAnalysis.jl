using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
import ElectroPhysiology: make_circular_roi!
using Plots
# using Pkg; Pkg.activate("test")
# using GLMakie, PhysiologyPlotting

#%% ╔═╡Found an issue. With really crazy spikes, the baseline correction is not working.
img_fn3 = raw"F:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b5_grabda-nircat-300uA_pulse014.tif"
stim_fn3 = raw"F:\Data\Patching\2025-05-15-GRAB-DA-STR\25515021.abf"

data_img = readImage(img_fn3)
deinterleave!(data_img)
stimulus = readABF(stim_fn3, stimulus_name = "IN 3", stimulus_threshold = 0.5, flatten_episodic = true)
spike_train_group!(stimulus, 1.0)
addStimulus!(data_img, stimulus, "IN 3")
truncate_data!(data_img, 0.0, 400.0)

#%% Our first 
stimulus_idx = 2
main_t_stim = getStimulusEndTime(data_img)[stimulus_idx]
trunc_start = main_t_stim - 40
trunc_end = main_t_stim + 120

stim_2 = truncate_data(data_img, trunc_start, trunc_end)
t_stim = getStimulusEndTime(stim_2)[1]
stim_frame = round(Int, t_stim ./ stim_2.dt)

pixel_splits_roi!(stim_2, 8)
roi_indices = getROImask(stim_2) |> unique
sig_rois = process_rois(stim_2, 
    stim_frame = stim_frame,
    window = 40,
    baseline_divisor_start = 20,
    baseline_divisor_end = 5,
    linear_fill_start = 5,
    linear_fill_end = 100,

    pos_sig_level = 2.0,
    neg_sig_level = 3.0,

    sig_threshold_std_start = 1,
    sig_threshold_std_end = 5,
    sig_threshold_mean_start = 12,
    sig_threshold_mean_end = 2,
    argmax_threshold_end = 100,
    max_dfof_end = 5,
    min_dfof_end = 100
)

stim_2.HeaderDict["sig_rois_idxs"]
stim_2.HeaderDict["sig_rois_mask_segment"]
sig_arr = getROIarr(stim_2, stim_2.HeaderDict["sig_rois_idxs"])
sig_arr_zproj = mean(sig_arr, dims = 1)[1,:,:]
sig_mask = reshape(sig_rois, sqrt(length(sig_rois)) |> Int64, sqrt(length(sig_rois)) |> Int64)
heatmap(rotl90(sig_mask), aspect_ratio = 1)

#%%
p1 = plot(sig_arr_zproj[:,1])
p2 = plot(sig_arr_zproj[:,2])

plot(p1, p2, layout = (2,1))

#%%
#Now lets draw a circular ROI and do the same calculations
make_circular_roi!(stim_2, 100, 100, 100)
