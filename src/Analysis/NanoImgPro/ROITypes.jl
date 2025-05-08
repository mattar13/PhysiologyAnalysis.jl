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
function get_significant_rois(analysis::ROIAnalysis, stim_idx::Int=1, channel_idx::Int=1)
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
    get_dfof_traces(analysis::ROIAnalysis, roi_id::Union{Int,Nothing}=nothing, stim_idx::Int=1, channel_idx::Int=1)

Get dF/F traces for specified ROI(s), stimulus, and channel.
"""
function get_dfof_traces(analysis::ROIAnalysis, roi_id::Union{Int,Nothing}=nothing, stim_idx::Int=1, channel_idx::Int=1)
    if isnothing(roi_id)
        return [filter(t -> t.channel == channel_idx, traces)[stim_idx].dfof for traces in values(analysis.rois)]
    end
    traces = filter(t -> t.channel == channel_idx, analysis.rois[roi_id])
    return traces[stim_idx].dfof
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