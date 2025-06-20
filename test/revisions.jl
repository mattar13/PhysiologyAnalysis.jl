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
println("Loading the quinpirole baseline data...")
main_channel = :red
img_fn3 = raw"G:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b5_grabda-nircat-300uA_pulse014.tif"
stim_fn3 = raw"G:\Data\Patching\2025-05-15-GRAB-DA-STR\25515021.abf"

data_img = readImage(img_fn3)
deinterleave!(data_img)
truncate_data!(data_img, t_begin = 0.0, t_end = 300.0)
stimulus = readABF(stim_fn3, stimulus_name = "IN 3", stimulus_threshold = 0.5, flatten_episodic = true)
spike_train_group!(stimulus, 1.0)
addStimulus!(data_img, stimulus, "IN 3")

function plot_baseline_trace(trace::AbstractVector{T}; 
    window::Int = 20,
    baseline_divisor_start = 1, baseline_divisor_end = nothing,
    linear_fill_start = nothing, linear_fill_end = nothing,
    kwargs...
) where T<:Real
    dFoF = baseline_trace(trace; window = window, baseline_divisor_start = baseline_divisor_start, baseline_divisor_end = baseline_divisor_end, linear_fill_start = linear_fill_start, linear_fill_end = linear_fill_end, kwargs...)
    p1 = plot(trace, label = "Raw")
    vline!(p1, [baseline_divisor_start, baseline_divisor_end], color = :red, label = "Baseline Divisor")
    vline!(p1, [linear_fill_start, linear_fill_end], color = :blue, label = "Linear Fill")
    p2 = plot(dFoF, label = "Filtered")
    plot(p1, p2, layout = (2, 1))
end

#Demonstrate the ROI significance finding

#%% Lets build the pixel splits function again (using minimal AI)
pixel_splits_roi!(data_img, 8)
mask = getROImask(data_img)
stim_end_times = getStimulusEndTime(data_img)
stim_end_idxs = round.(Int, stim_end_times ./ data_img.dt)

sig_rois, dFoF_traces = process_rois(data_img, stim_frame = stim_end_idxs[2])
t_rng = data_img.t[1]:data_img.dt:data_img.t[end]
roi_trace = mean(getROIarr(data_img, 1), dims = 1)[1,:, 2]

sig_mask = reshape(sig_rois, sqrt(length(sig_rois)) |> Int64, sqrt(length(sig_rois)) |> Int64)
heatmap(sig_mask)

window::Int=40
baseline_divisor_start = 2
baseline_divisor_end = nothing
linear_fill_start = nothing
linear_fill_end = nothing

dFoF = baseline_trace(roi_trace, 
    window = window, 
    baseline_divisor_start = baseline_divisor_start, 
    baseline_divisor_end = baseline_divisor_end, 
    linear_fill_start = linear_fill_start, 
    linear_fill_end = linear_fill_end
)


#%% Plot everything
pos_sig_level = 2.0
neg_sig_level = 3.0
stim_frame = stim_end_idxs[2]
sig_threshold_std_start = 1
sig_threshold_std_end = 5
sig_threshold_mean_start = 12
sig_threshold_mean_end = 2
argmax_threshold_end = 25
max_dfof_end = 5
min_dfof_end = 100

#These are the parameters for the signal threshold
sig_std = std(dFoF[sig_threshold_std_start:stim_frame-sig_threshold_std_end])
sig_mean = mean(dFoF[stim_frame-sig_threshold_mean_start: stim_frame-sig_threshold_mean_end])
#println("sig_mean: $sig_mean, sig_std: $sig_std")

#Calculate the threshold for the signal
sig_threshold = sig_mean + sig_std * pos_sig_level
neg_threshold = sig_mean - sig_std * neg_sig_level

trace_argmax = argmax(dFoF[stim_frame:stim_frame+argmax_threshold_end]) + stim_frame
        
max_dfof = mean(dFoF[trace_argmax:trace_argmax+max_dfof_end])
min_dfof = minimum(dFoF[stim_frame:stim_frame+min_dfof_end])
over_max = max_dfof > sig_threshold
over_min = min_dfof > neg_threshold

p1 = plot(t_rng, roi_trace, label = "Raw")
p2 = plot(t_rng, dFoF, label = "Filtered")
vline!(p2, getStimulusEndTime(data_img), color = :red, label = "Stimulus")
plot(p1, p2, layout = (2, 1))

#%%

