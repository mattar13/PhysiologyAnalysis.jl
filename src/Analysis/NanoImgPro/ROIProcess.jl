"""
Perform overall baseline correction using ALS and centered Moving Average.

- `stim_frame=0`: Disables stimulus-based normalization.
- `ma_tail_start=0`: Disables ignored region in the moving average correction.
- `window`: Size of the moving average window.

Returns a baseline-corrected dF/F trace.
"""
function baseline_trace(trace::AbstractVector{T}; 
    stim_frame = nothing, window::Int=15, #This is the window of the moving average for dF
    kwargs...
) where T<:Real

    # Normalize using pre-stimulus baseline if stim_frame is set
    if isnothing(stim_frame)
        baseline_divisor = mean(trace)  # Default to global mean if no stimulus
    else
        # Calculate the mean of the trace before the stimulus frame
        baseline_divisor = mean(trace[1:stim_frame])
    end
    F0 = trace ./ baseline_divisor

    # Apply Asymmetric Least Squares (ALS) smoothing
    drift = baseline_als(F0; kwargs...)
    dF = F0 .- drift

    dFoF = moving_average(dF; window = window)

    return dFoF
end

#Baseline trace should take in only the trace and window size, but when we add stims, we can properly Adjust
function roi_processing(data::Experiment{TWO_PHOTON, T}, roi_index, stim_index; 
    channel = nothing, delay_time = 20.0, 
    stim_frame = nothing, window::Int=15, #This is the window of the moving average for dF
    verbose = false
    kwargs...
) where T<:Real

    data_dict = Dict{Symbol, Any}()
    stims = getStimulusEndTime(data)
    if isnothing(channel)
        for i in eachchannel(data)
            println("Processing channel $i")
            channel_dict = roi_processing(data, roi_index, stim_index; channel = i, stim_frame = stim_frame, window = window, kwargs...)
            if isempty(data_dict)
                data_dict = channel_dict
            else
                data_dict["roi_trace"] = vcat(data_dict["roi_trace"], channel_dict["roi_trace"])
                data_dict["dFoF"] = vcat(data_dict["dFoF"], channel_dict["dFoF"])
            end
        end
        return data_dict
    else    
        if verbose
            println("Processing stimulus indices for channel $channel")
        end
        data_dict["stim_start_time"] = stim_start_time = getStimulusEndTime(data)[stim_index] - delay_time
        data_dict["stim_end_time"] = stim_end_time = getStimulusStartTime(data)[stim_index+1]
        data_dict["stim_start_index"] = stim_start_index = round(Int64, (stim_start_time)/data.dt)
        data_dict["stim_end_index"] = stim_end_index = round(Int64, (stim_end_time)/data.dt)
        
        roi_frames = getROIarr(data, roi_index) #This is the ROI array for the first movie
        data_dict["roi_trace"] = roi_trace = mean(roi_frames, dims = (1))[1, stim_start_index:stim_end_index, channel] #This is the mean of the ROI array
        
        data_dict["dFoF"] = dFoF = baseline_trace(roi_trace; stim_frame = stim_frame, window = window, kwargs...)
        
        return data_dict
    end
    #Next we want to baseline the trace
end