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
"""

"""
function roi_processing(data::Experiment{TWO_PHOTON, T}, stim_index; 
    roi_index = nothing, channel_index = nothing, 
    delay_time = 150.0, window::Int=15, #This is the window of the moving average for dF
    verbose = false, n_stds = 2.0, #This is the number of standard deviations to use for the thresholding
    kwargs...
) where T<:Real

    data_dict = Dict{String, Any}()
    stims = getStimulusEndTime(data)
    roi_list = getROImask(data) |> unique |> sort

    if isnothing(channel_index) && isnothing(roi_index)
        for roi_idx in roi_list
            println("Processing roi: $roi_idx")
            channel_dict = roi_processing(data, stim_index; roi_index = roi_idx, 
                window = window, n_stds, kwargs...
            )
            if isempty(data_dict)
                data_dict = [channel_dict]
            else
                push!(data_dict, channel_dict)
            end
        end
        return data_dict
    elseif isnothing(channel_index)
        #We should deal with this too
        for ch_idx in axes(data, 3)
            #println("Processing channel $ch_idx")
            
            channel_dict = roi_processing(data, stim_index; 
                roi_index = roi_index, channel_index = ch_idx, 
                window = window, n_stds = n_stds, kwargs...
            )

            if isempty(data_dict)
                data_dict = channel_dict
            else
                data_dict["roi_trace"] = hcat(data_dict["roi_trace"], channel_dict["roi_trace"])
                data_dict["dFoF"] = hcat(data_dict["dFoF"], channel_dict["dFoF"])
                data_dict["fit_param"] = hcat(data_dict["fit_param"], channel_dict["fit_param"])
                data_dict["significance"] = hcat(data_dict["significance"], channel_dict["significance"])
            end
        end
        return data_dict
    elseif isnothing(roi_index)
        for roi_idx in roi_list
            println("Processing roi: $roi_idx")
            channel_dict = roi_processing(data, stim_index; roi_index = roi_idx, 
                channel_index = channel_index,
                window = window, n_stds, kwargs...
            )
            if isempty(data_dict)
                data_dict = [channel_dict]
            else
                push!(data_dict, channel_dict)
            end
        end
        return data_dict
    else
        if verbose
            println("Processing stimulus indices for channel $channel")
        end
        data_dict["stim_start_time"] = stim_start_time = max(getStimulusEndTime(data)[stim_index] - delay_time, 0.0)
        data_dict["stim_end_time"] = stim_end_time = min(getStimulusStartTime(data)[stim_index+1], data.t[end])

        data_dict["stim_start_index"] = stim_start_index = round(Int64, (stim_start_time+1)/data.dt)
        data_dict["stim_end_index"] = stim_end_index = round(Int64, (stim_end_time)/data.dt)
        
        roi_frames = getROIarr(data, roi_index) #This is the ROI array for the first movie
        data_dict["roi_trace"] = roi_trace = mean(roi_frames, dims = (1))[1, stim_start_index:stim_end_index, channel_index] #This is the mean of the ROI array
        
        data_dict["dFoF"] = dFoF = baseline_trace(roi_trace; stim_frame = stim_start_index, window = window, kwargs...)
        data_dict["t_series"] = t_series = collect(1:length(dFoF))*data.dt
        
        #This is similar to OTSU thresholding
        sig_region = dFoF[1:stim_end_index] #We can choose
        pos_sig_threshold = mean(sig_region) + n_stds*std(sig_region)
        # neg_sig_threshold = mean(sig_region) - n_stds*std(sig_region)

        #Check if there are any values above this threshold after the stim_end_index
        max_post_stim_res = maximum(dFoF[stim_start_index:stim_end_index])
        if max_post_stim_res < pos_sig_threshold
            data_dict["significance"] = false
        else
            data_dict["significance"] = true
        end

        #Fit the data to the parametric model
        lb = [max(maximum(dFoF)-0.5, 0.0), 0.0, 0.0, max(delay_time-10.0, 0.0)] #The shift time should be the last param
        p0 = [maximum(dFoF), 100.0, 10.0, delay_time] #This is the initial guess for the fit
        ub = [maximum(dFoF)+1.0, 1000.0, 1000.0, delay_time+10.0] #This is the upper bound for the fit
        fit = fit_parametric(t_series, dFoF, lb = lb, p0 = p0, ub = ub)
        data_dict["fit_param"] = fit.param

        return data_dict
    end
end

#Now write some easy extraction functions to get the data out of the dictionary
function get_significance(data_dict::Array{Dict{String, Any}})
    map(dd->dd["significance"], data_dict)
end

function get_dFoF(data_dict::Array{Dict{String, Any}})
    vcat(map(dd->dd["dFoF"], data_dict)'...)
end

function get_fit_param(data_dict::Array{Dict{String, Any}})
    map(dd->dd["fit_param"], data_dict)
end