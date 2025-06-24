using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
import ElectroPhysiology: make_circular_roi!
using GLMakie
GLMakie.activate!()
using JSON
# using Pkg; Pkg.activate("test")
# using GLMakie, PhysiologyPlotting

#===============================================#
#%%Load parameters from JSON file
#===============================================#

# Load parameters from the JSON file
json_path = raw"C:\Users\Matt\JuliaDev\ElectroPhysiologyWorkspace\DataAnalysis\SWCNT_Datasets\Manuscript Figures\parameters.json"
json_string = read(json_path, String)
parameters = JSON.parse(json_string)

# Extract parameters by category
trunc_params = parameters["truncation"]
baseline_params = parameters["baselining"]
roi_params = parameters["roi_processing"]

# Assign parameters to variables for backward compatibility
stimulus_idx = trunc_params["stimulus_idx"]
trunc_before_stim = trunc_params["trunc_before_stim"]
trunc_after_stim = trunc_params["trunc_after_stim"]

window = baseline_params["window"]
baseline_divisor_start = baseline_params["baseline_divisor_start"]
baseline_divisor_end = baseline_params["baseline_divisor_end"]
linear_fill_start = baseline_params["linear_fill_start"]
linear_fill_end = baseline_params["linear_fill_end"]

pos_sig_level = roi_params["pos_sig_level"]
neg_sig_level = roi_params["neg_sig_level"]

sig_threshold_std_start = roi_params["sig_threshold_std_start"]
sig_threshold_std_end = roi_params["sig_threshold_std_end"]
sig_threshold_mean_start = roi_params["sig_threshold_mean_start"]
sig_threshold_mean_end = roi_params["sig_threshold_mean_end"]
argmax_threshold_end = roi_params["argmax_threshold_end"]
max_dfof_end = roi_params["max_dfof_end"]
min_dfof_end = roi_params["min_dfof_end"]

#Extract all of the stimulus indices


#===============================================#
#Load the data
#===============================================#
drive_root = "G:"
# Filenames for the analysis, using the b5 dataset from openingData.jl
img_fn = "$(drive_root)\\Data\\Two Photon\\2025-05-15-GRAB-DA_STR\\b5_grabda-nircat-300uA_pulse014.tif"
stim_fn = "$(drive_root)\\Data\\Patching\\2025-05-15-GRAB-DA-STR\\25515021.abf"

data_img = readImage(img_fn)
deinterleave!(data_img)
stimulus = readABF(stim_fn, stimulus_name = "IN 3", stimulus_threshold = 0.5, flatten_episodic = true)
spike_train_group!(stimulus, 1.0)
addStimulus!(data_img, stimulus, "IN 3")
truncate_data!(data_img, 0.0, 300.0)

main_t_stim = getStimulusEndTime(data_img)[stimulus_idx]
trunc_start = main_t_stim - trunc_before_stim
trunc_end = main_t_stim + trunc_after_stim

# --- GLMakie Figure Setup --- #
fig = Figure(size = (800, 800))

# --- Plot 1: Raw Trace with Stimulus and Truncation --- #
ax0 = Axis(fig[1, 1], title = "Raw Trace", xlabel = "Time (s)", ylabel = "Mean Intensity")
lines!(ax0, data_img.t, mean(data_img.data_array[:,:,2], dims = 1)[1,:], color = :black)
vlines!(ax0, getStimulusEndTime(data_img), color = :red, label = "Stimulus", linewidth = 2, alpha = 0.2)
vlines!(ax0, [main_t_stim], color = :green, label = "Main Stimulus", linewidth = 2)
vspan!(ax0, trunc_start, trunc_end, color = (:blue, 0.2))
axislegend(ax0)

#===============================================#
#Conduct the baseline correction
#===============================================#
pixel_ROI_data = truncate_data(data_img, trunc_start, trunc_end)
circular_ROI_data = copy(pixel_ROI_data)
t_stim = getStimulusEndTime(pixel_ROI_data)[1]
current_stim_frame = stim_frame = round(Int, t_stim ./ pixel_ROI_data.dt)
ir_img = mean(pixel_ROI_data.data_array[:,:,2], dims = 1)[1,:]

println("Current Stim Frame: $current_stim_frame")

linear_fill_start_idx = max(1, stim_frame - linear_fill_start)
linear_fill_end_idx = min(size(pixel_ROI_data, 2), stim_frame + linear_fill_end)
baseline_divisor_start_idx = max(1, stim_frame - baseline_divisor_start)
baseline_divisor_end_idx = max(1, stim_frame - baseline_divisor_end)

# --- Plot 2: Baseline Correction --- #
ax1 = Axis(fig[2, 1], title = "Baseline Correction", xlabel = "Time (s)", ylabel = "Mean Intensity")
lines!(ax1, pixel_ROI_data.t, ir_img, color = :black)
vlines!(ax1, getStimulusEndTime(pixel_ROI_data), color = :green, label = "Stimulus", linewidth = 2)
# Baseline divisor section
t_baseline_start = pixel_ROI_data.t[baseline_divisor_start_idx]
t_baseline_end = pixel_ROI_data.t[baseline_divisor_end_idx]
vspan!(ax1, t_baseline_start, t_baseline_end, color = (:red, 0.2), label = "Baseline Divisor")
# Linear fill section
t_linear_fill_start = pixel_ROI_data.t[linear_fill_start_idx]
t_linear_fill_end = pixel_ROI_data.t[linear_fill_end_idx]
vspan!(ax1, t_linear_fill_start, t_linear_fill_end, color = (:blue, 0.2), label = "Linear Fill")
axislegend(ax1, postion = :rb)

dFoF = baseline_trace(ir_img;
    window = window,
    stim_frame = stim_frame,
    baseline_divisor_start = baseline_divisor_start,
    baseline_divisor_end = baseline_divisor_end,
    linear_fill_start = linear_fill_start,
    linear_fill_end = linear_fill_end
)

#===============================================#
#Find the signal threshold
#===============================================#
println("Current Stim Frame: $current_stim_frame")
println(sig_threshold_std_start)

data = pixel_ROI_data
sig_threshold_std_start_idx = sig_threshold_std_start
sig_threshold_std_end_idx = min(size(data, 2), current_stim_frame - sig_threshold_std_end)
sig_threshold_mean_start_idx = max(1, current_stim_frame - sig_threshold_mean_start)
sig_threshold_mean_end_idx = min(size(data, 2), current_stim_frame - sig_threshold_mean_end)

sig_std = std(dFoF[sig_threshold_std_start_idx:sig_threshold_std_end_idx])
sig_mean = mean(dFoF[sig_threshold_mean_start_idx: sig_threshold_mean_end_idx])
sig_threshold = sig_mean + sig_std * pos_sig_level
neg_threshold = sig_mean - sig_std * neg_sig_level
# --- Plot 3: dF/F and Signal Thresholds --- #
ax2 = Axis(fig[3, 1], title = "dF/F and Signal Thresholds", xlabel = "Time (s)", ylabel = "dF/F")
lines!(ax2, pixel_ROI_data.t, dFoF, color = :black)
vlines!(ax2, [t_stim], color = :green, label = "Stimulus", linewidth = 2)
hlines!(ax2, [sig_threshold, neg_threshold], label = "Signal Threshold", linestyle = :dash, color = :blue, linewidth = 2, alpha = 0.7)

argmax_threshold_end_idx = min(size(data, 2), current_stim_frame + argmax_threshold_end)
max_dfof_end_idx = max(1, current_stim_frame + max_dfof_end)
min_dfof_end_idx = min(size(data, 2), current_stim_frame + min_dfof_end)

t_min_dfof_start = pixel_ROI_data.t[current_stim_frame]
t_min_dfof_end = pixel_ROI_data.t[min_dfof_end_idx]
trace_argmax = argmax(dFoF[current_stim_frame:argmax_threshold_end_idx]) + current_stim_frame
t_max_dfof_start = pixel_ROI_data.t[trace_argmax]
t_max_dfof_end = pixel_ROI_data.t[max_dfof_end_idx]

max_dfof = mean(dFoF[trace_argmax:trace_argmax+max_dfof_end])
min_dfof = minimum(dFoF[current_stim_frame:min_dfof_end_idx])

hlines!(ax2, [max_dfof, min_dfof], color = :orange, label = "Max/Min DF/F", linewidth = 2, alpha = 0.5)
vspan!(ax2, t_min_dfof_start, t_min_dfof_end, color = (:blue, 0.2), label = "Min DFF region")
vspan!(ax2, t_max_dfof_start, t_max_dfof_end, color = (:green, 0.2), label = "Max DFF region")
axislegend(ax2)

fig
#===============================================#
#%%Find the ROIs
#===============================================#
nx, ny = getIMG_size(pixel_ROI_data)
pixel_splits_roi!(pixel_ROI_data, 8)
roi_indices = getROImask(pixel_ROI_data) |> unique
sig_rois = process_rois(pixel_ROI_data, 
    pos_sig_level = pos_sig_level,
    neg_sig_level = neg_sig_level,
    
    window = window,
    
    baseline_divisor_start = baseline_divisor_start,
    baseline_divisor_end = baseline_divisor_end,
    linear_fill_start = linear_fill_start,
    linear_fill_end = linear_fill_end,

    sig_threshold_std_start = sig_threshold_std_start,
    sig_threshold_std_end = sig_threshold_std_end,
    sig_threshold_mean_start = sig_threshold_mean_start,
    sig_threshold_mean_end = sig_threshold_mean_end,
    argmax_threshold_end = argmax_threshold_end,
    max_dfof_end = max_dfof_end,
    min_dfof_end = min_dfof_end
)

sig_arr = getROIarr(pixel_ROI_data, pixel_ROI_data.HeaderDict["sig_rois_idxs"])
sig_arr_zproj = mean(sig_arr, dims = 1)[1,:,:]
sig_mask = reshape(sig_rois, nx÷8, ny÷8)
sig_arr_dFoF = zeros(size(sig_arr_zproj))

sig_arr_dFoF[:,1] = baseline_trace(sig_arr_zproj[:,1],
    window = window,
    stim_frame = stim_frame,
    baseline_divisor_start = baseline_divisor_start,
    baseline_divisor_end = baseline_divisor_end,
    linear_fill_start = linear_fill_start,
    linear_fill_end = linear_fill_end
)
sig_arr_dFoF[:,2] = baseline_trace(sig_arr_zproj[:,2],
    window = window,
    stim_frame = stim_frame,
    baseline_divisor_start = baseline_divisor_start,
    baseline_divisor_end = baseline_divisor_end,
    linear_fill_start = linear_fill_start,
    linear_fill_end = linear_fill_end
)
# --- Plot 4: ROI Heatmap and Traces --- #
ax3a = Axis(fig[1,2], title = "Significant ROI Mask", aspect = 1)
heatmap!(ax3a, rotl90(sig_mask), colormap = :viridis)
ax3b = Axis(fig[2,2], title = "ROI Trace 1", xlabel = "Pixel Index", ylabel = "Mean Intensity")
lines!(ax3b, sig_arr_dFoF[:,1], color = :blue)
ax3c = Axis(fig[3,2], title = "ROI Trace 2", xlabel = "Pixel Index", ylabel = "Mean Intensity")
lines!(ax3c, sig_arr_dFoF[:,2], color = :red)

fig