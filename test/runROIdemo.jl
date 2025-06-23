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

#===============================================#
#%%Define the parameters for the ROI significance finding
#===============================================#

stimulus_idx = 2
trunc_before_stim = 50
trunc_after_stim = 120

window = 40
baseline_divisor_start = 20
baseline_divisor_end = 0
linear_fill_start = 5
linear_fill_end = 100

pos_sig_level = 2.0
neg_sig_level = 3.0

sig_threshold_std_start = 1
sig_threshold_std_end = 5
sig_threshold_mean_start = 12
sig_threshold_mean_end = 2
argmax_threshold_end = 25
max_dfof_end = 25
min_dfof_end = 100

#===============================================#
#%%Load the data
#===============================================#
img_fn3 = raw"F:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b5_grabda-nircat-300uA_pulse014.tif"
stim_fn3 = raw"F:\Data\Patching\2025-05-15-GRAB-DA-STR\25515021.abf"

data_img = readImage(img_fn3)
deinterleave!(data_img)
stimulus = readABF(stim_fn3, stimulus_name = "IN 3", stimulus_threshold = 0.5, flatten_episodic = true)
spike_train_group!(stimulus, 1.0)
addStimulus!(data_img, stimulus, "IN 3")
truncate_data!(data_img, 0.0, 300.0)

main_t_stim = getStimulusEndTime(data_img)[stimulus_idx]
trunc_start = main_t_stim - trunc_before_stim
trunc_end = main_t_stim + trunc_after_stim

p0 = plot(data_img.t, mean(data_img.data_array[:,:,2], dims = 1)[1,:])
vline!(p0, getStimulusEndTime(data_img), color = :red, label = "Stimulus", alpha = 0.2)
vline!(p0, [main_t_stim], color = :green, label = "Main Stimulus")
vspan!(p0, [trunc_start, trunc_end], color = :blue, label = "Truncation Range", alpha = 0.2)

#===============================================#
#%%Conduct the baseline correction
#===============================================#
pixel_ROI_data = truncate_data(data_img, trunc_start, trunc_end)
circular_ROI_data = copy(pixel_ROI_data)
t_stim = getStimulusEndTime(pixel_ROI_data)[1]
stim_frame = round(Int, t_stim ./ pixel_ROI_data.dt)
ir_img = mean(pixel_ROI_data.data_array[:,:,2], dims = 1)[1,:]

p1 = plot(pixel_ROI_data.t, ir_img, legend = :topright)
vline!(p1, getStimulusEndTime(pixel_ROI_data), color = :green, label = "Stimulus")

#Plot the baseline divisor section
t_baseline_start = pixel_ROI_data.t[stim_frame - baseline_divisor_start]
t_baseline_end = pixel_ROI_data.t[stim_frame - baseline_divisor_end]
vspan!(p1, [t_baseline_start, t_baseline_end], color = :red, label = "Baseline Divisor", alpha = 0.2)

#Plot the linear fill section
t_linear_fill_start = pixel_ROI_data.t[stim_frame - linear_fill_start]
t_linear_fill_end = pixel_ROI_data.t[stim_frame + linear_fill_end]
vspan!(p1, [t_linear_fill_start, t_linear_fill_end], color = :blue, label = "Linear Fill", alpha = 0.2)

dFoF = baseline_trace(ir_img;
    window = 40,
    stim_frame = stim_frame,
    baseline_divisor_start = baseline_divisor_start,
    baseline_divisor_end = baseline_divisor_end,
    linear_fill_start = linear_fill_start,
    linear_fill_end = linear_fill_end
)
p1

#===============================================#
#%%Find the signal threshold
#===============================================#

p2 = plot(pixel_ROI_data.t, dFoF)
vline!(p2, [t_stim], color = :green, label = "Stimulus")
sig_std = std(dFoF[sig_threshold_std_start:stim_frame-sig_threshold_std_end])
sig_mean = mean(dFoF[stim_frame-sig_threshold_mean_start: stim_frame-sig_threshold_mean_end])
sig_threshold = sig_mean + sig_std * pos_sig_level
neg_threshold = sig_mean - sig_std * neg_sig_level
hline!(p2, [sig_threshold, neg_threshold], color = :red, label = "Signal Threshold", alpha = 0.2)

#
t_min_dfof_start = pixel_ROI_data.t[stim_frame]
t_min_dfof_end = pixel_ROI_data.t[stim_frame + min_dfof_end]
trace_argmax = argmax(dFoF[stim_frame:stim_frame+argmax_threshold_end]) + stim_frame
t_max_dfof_start = pixel_ROI_data.t[trace_argmax]
t_max_dfof_end = pixel_ROI_data.t[trace_argmax + max_dfof_end]

max_dfof = mean(dFoF[trace_argmax:trace_argmax+max_dfof_end])
min_dfof = minimum(dFoF[stim_frame:stim_frame+min_dfof_end])
hline!(p2, [max_dfof, min_dfof], color = :orange, label = "Max/Min DF/F", alpha = 0.5)

vspan!(p2, [t_min_dfof_start, t_min_dfof_end], color = :blue, label = "Min DF/F", alpha = 0.2)
vspan!(p2, [t_max_dfof_start, t_max_dfof_end], color = :green, label = "Max DF/F", alpha = 0.2)

plt_sig_find = plot(p0, p1, p2, layout = (3, 1))

#===============================================#
#%%Find the ROIs
#===============================================#

pixel_splits_roi!(pixel_ROI_data, 8)
roi_indices = getROImask(pixel_ROI_data) |> unique
sig_rois = process_rois(pixel_ROI_data, stim_idx =   1)

pixel_ROI_data.HeaderDict["sig_rois_idxs"]
pixel_ROI_data.HeaderDict["sig_rois_mask_segment"]
sig_arr = getROIarr(pixel_ROI_data, pixel_ROI_data.HeaderDict["sig_rois_idxs"])
sig_arr_zproj = mean(sig_arr, dims = 1)[1,:,:]
sig_mask = reshape(sig_rois, sqrt(length(sig_rois)) |> Int64, sqrt(length(sig_rois)) |> Int64)

hm = heatmap(rotl90(sig_mask), aspect_ratio = 1, c = :viridis)
p1 = plot(sig_arr_zproj[:,1])
p2 = plot(sig_arr_zproj[:,2])

plot(hm, hm, p1, p2, layout = (2,2))