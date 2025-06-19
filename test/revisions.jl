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
truncate_data!(data_img, t_begin = 0.0, t_end = 500.0)
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
mask_xy = round.(Int64, sqrt(size(sig_rois,1))) 
heatmap(reshape(sig_rois, mask_xy, mask_xy))

#%% We want to store the ROIs in a object
#make_circular_roi!(data_img, (75, 75), 75)

t_rng = data_img.t[1]:data_img.dt:data_img.t[end]
roi_trace = mean(getROIarr(data_img, 50), dims = 1)[1,:, 2]

roi_FILT = baseline_trace(roi_trace, 
    window = 40, 
    baseline_divisor_start = 1, baseline_divisor_end = 120,
    linear_fill_start = 120, linear_fill_end = 180
)

p1 = plot(t_rng, roi_trace, label = "Raw")
p2 = plot(t_rng, roi_FILT, label = "Filtered")
vline!(p2, getStimulusEndTime(data_img), color = :red, label = "Stimulus")
plot(p1, p2, layout = (2, 1))

#%%

