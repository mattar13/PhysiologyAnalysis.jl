"""
Types for representing ROI analysis results from two-photon imaging data.
"""

"""
    ROITrace

Represents the analysis results for a single ROI trace.

Fields:
- `raw_trace`: The raw fluorescence trace
- `dfof`: The ΔF/F trace after baseline correction
- `t_series`: Time points for the trace
- `stim_start_time`: Start time of the stimulus (s)
- `stim_end_time`: End time of the stimulus (s)
- `channel`: Channel index
- `stimulus_index`: Index of the stimulus this trace corresponds to
- `fit_parameters`: Parameters from curve fitting [amplitude, tau_on, tau_off, delay]
- `is_significant`: Whether the response is significant
"""
struct ROITrace
    raw_trace::Vector{Float64}
    dfof::Vector{Float64}
    t_series::Vector{Float64}
    stim_start_time::Float64
    stim_end_time::Float64
    channel::Int
    stimulus_index::Int
    fit_parameters::Vector{Float64}
    is_significant::Bool
end

"""
    ROIAnalysis

Contains all ROI analysis results for a two-photon imaging experiment.

Fields:
- `rois`: Dictionary mapping ROI ID to vector of ROITrace objects (one per stimulus per channel)
- `stimulus_indices`: Vector of analyzed stimulus indices
- `channels`: Vector of analyzed channel indices
- `analysis_parameters`: Dictionary of analysis parameters (window size, delay time, etc.)
"""
struct ROIAnalysis
    rois::Dict{Int, Vector{ROITrace}}
    stimulus_indices::Vector{Int}
    channels::Vector{Int}
    analysis_parameters::Dict{Symbol, Any}
end

# Constructor with default analysis parameters
function ROIAnalysis(rois::Dict{Int, Vector{ROITrace}}, 
    stimulus_indices::Vector{Int}, 
    channels::Vector{Int};
    delay_time::Float64=50.0,
    window::Int=15,
    n_stds::Float64=2.0,
    kwargs...
)
    analysis_parameters = Dict{Symbol, Any}(
        :delay_time => delay_time,
        :window => window,
        :n_stds => n_stds
    )
    # Add any additional parameters from kwargs
    merge!(analysis_parameters, Dict(Symbol(k)=>v for (k,v) in kwargs))
    
    ROIAnalysis(rois, stimulus_indices, channels, analysis_parameters)
end

# Utility functions
"""
    get_significant_rois(analysis::ROIAnalysis, stim_idx::Int=1, channel_idx::Int=1)

Get indices of ROIs with significant responses for a given stimulus and channel.
"""
function get_significant_rois(analysis::ROIAnalysis; channel_idx::Int=1)
    return [id for (id, traces) in analysis.rois if any(t -> t.channel == channel_idx && t.is_significant, traces)]
end

"""
    get_mean_response(analysis::ROIAnalysis, roi_id::Int, channel_idx::Int=1)

Calculate the mean dF/F response across all stimuli for a given ROI and channel.
"""
function get_mean_response(analysis::ROIAnalysis, roi_id::Int, channel_idx::Int=1)
    traces = filter(t -> t.channel == channel_idx, analysis.rois[roi_id])
    return mean(t.dfof for t in traces)
end

"""
    get_mean_response(analysis::ROIAnalysis, roi_ids::Vector{Int}, channel_idx::Int=1)

Calculate the mean dF/F response across all stimuli for multiple ROIs and channel.
"""
function get_mean_response(analysis::ROIAnalysis, roi_ids::Vector{Int}, channel_idx::Int=1)
    return [get_mean_response(analysis, roi_id, channel_idx) for roi_id in roi_ids]
end

"""
    get_roi_traces(analysis::ROIAnalysis, roi_id::Union{Int,Nothing}=nothing, stim_idx::Int=1, channel_idx::Int=1)

Get raw traces for specified ROI(s), stimulus, and channel.
"""
function get_roi_traces(analysis::ROIAnalysis, roi_id::Union{Int,Nothing}=nothing, stim_idx::Int=1, channel_idx::Int=1)
    if isnothing(roi_id)
        return [filter(t -> t.channel == channel_idx, traces)[stim_idx].raw_trace for traces in values(analysis.rois)]
    end
    traces = filter(t -> t.channel == channel_idx, analysis.rois[roi_id])
    return traces[stim_idx].raw_trace
end

"""
    get_roi_traces(analysis::ROIAnalysis, roi_ids::Vector{Int}, stim_idx::Int=1, channel_idx::Int=1)

Get raw traces for multiple ROIs, stimulus, and channel.
"""
function get_roi_traces(analysis::ROIAnalysis, roi_ids::Vector{Int}, stim_idx::Int=1, channel_idx::Int=1)
    return [get_roi_traces(analysis, roi_id, stim_idx, channel_idx) for roi_id in roi_ids]
end

"""
    get_dfof_traces(analysis::ROIAnalysis, roi_id::Union{Int,Nothing}=nothing, stim_idx::Int=1, channel_idx::Int=1)

Get dF/F traces for specified ROI(s), stimulus, and channel.
"""
function get_dfof_traces(analysis::ROIAnalysis, roi_id::Union{Int,Nothing}=nothing; stim_idx::Int=1, channel_idx::Int=1)
    if isnothing(roi_id)
        return [filter(t -> t.channel == channel_idx, traces)[stim_idx].dfof for traces in values(analysis.rois)]
    end
    traces = filter(t -> t.channel == channel_idx, analysis.rois[roi_id])
    return traces[stim_idx].dfof
end

"""
    get_dfof_traces(analysis::ROIAnalysis, roi_ids::Vector{Int}, stim_idx::Int=1, channel_idx::Int=1)

Get dF/F traces for multiple ROIs, stimulus, and channel.
"""
function get_dfof_traces(analysis::ROIAnalysis, roi_ids::Vector{Int}; stim_idx::Int=1, channel_idx::Int=1)
    return [get_dfof_traces(analysis, roi_id; stim_idx = stim_idx, channel_idx = channel_idx) for roi_id in roi_ids]
end

"""
    get_fit_parameters(analysis::ROIAnalysis, roi_id::Union{Int,Nothing}=nothing, stim_idx::Int=1, channel_idx::Int=1)

Get fit parameters for specified ROI(s), stimulus, and channel.
"""
function get_fit_parameters(analysis::ROIAnalysis, roi_id::Union{Int,Nothing}=nothing, stim_idx::Int=1, channel_idx::Int=1)
    if isnothing(roi_id)
        return [filter(t -> t.channel == channel_idx, traces)[stim_idx].fit_parameters for traces in values(analysis.rois)]
    end
    traces = filter(t -> t.channel == channel_idx, analysis.rois[roi_id])
    return traces[stim_idx].fit_parameters
end

"""
    get_fit_parameters(analysis::ROIAnalysis, roi_ids::Vector{Int}, stim_idx::Int=1, channel_idx::Int=1)

Get fit parameters for multiple ROIs, stimulus, and channel.
"""
function get_fit_parameters(analysis::ROIAnalysis, roi_ids::Vector{Int}, stim_idx::Int=1, channel_idx::Int=1)
    return [get_fit_parameters(analysis, roi_id, stim_idx, channel_idx) for roi_id in roi_ids]
end

"""
    get_significant_traces_by_channel(analysis)

Return a vector of matrices, one per channel, each of size (n_ROIs, n_Datapoints).
All traces are truncated to the global minimum length for consistency.
"""
function get_significant_traces_by_channel(analysis)
    n_channels = length(analysis.channels)
    traces_by_channel = Vector{Matrix{Float64}}(undef, n_channels)
    min_length = typemax(Int)
    # First, collect traces for each channel
    traces_per_channel = [Vector{Vector{Float64}}() for _ in 1:n_channels]
    for (_, traces) in analysis.rois
        for trace in traces
            if trace.is_significant
                push!(traces_per_channel[trace.channel], trace.dfof)
                min_length = min(min_length, length(trace.dfof))
            end
        end
    end
    # Now, build the matrices
    for ch in 1:n_channels
        n_rois = length(traces_per_channel[ch])
        if n_rois == 0 || min_length == 0
            traces_by_channel[ch] = zeros(0, 0)
        else
            mat = zeros(n_rois, min_length)
            for i in 1:n_rois
                mat[i, :] = traces_per_channel[ch][i][1:min_length]
            end
            traces_by_channel[ch] = mat
        end
    end
    return traces_by_channel
end

"""
    get_significant_roi_pixels(analysis::ROIAnalysis, data::Experiment{TWO_PHOTON}; channel_idx::Int=1)

Get the pixel coordinates for all significant ROIs in the original image.

Parameters:
- `analysis`: The ROIAnalysis object containing processed ROI data
- `data`: The original Experiment object containing the ROI mask
- `channel_idx`: Channel index to use for determining significant ROIs (default: 1)

Returns:
- Vector of linear indices for pixels in significant ROIs from the specified channel
"""
function get_significant_roi_pixels(analysis::ROIAnalysis, data::Experiment{TWO_PHOTON}; channel_idx::Int=1)
    # Get significant ROIs for the specified channel
    significant_rois = get_significant_rois(analysis, 1, channel_idx)
    
    # Get ROI mask from experiment
    roi_mask = getROImask(data)
    
    # Find all pixels belonging to any significant ROI
    pixel_coords = findall(roi_mask .∈ (significant_rois,))
    
    # Convert Cartesian indices to linear indices
    return LinearIndices(roi_mask)[pixel_coords]
end 