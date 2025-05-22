using ProgressMeter

"""
    process_rois(data::Experiment{TWO_PHOTON, T}; 
        stim_indices=nothing,
        channels=nothing,
        roi_indices=nothing,
        delay_time=50.0,
        window::Int=15,
        n_stds=2.0,
        sig_window=50.0,  # Time window in ms to look for significant responses after stimulus
        kwargs...
    ) where T<:Real

Process regions of interest (ROIs) from two-photon imaging data with the following steps:
1. Split pixels into ROIs if not already done (use pixel_splits_roi!)
2. For each ROI and channel:
   - Extract ROI traces
   - Calculate baseline using pre-stimulus period
   - Compute dF/F using moving average correction
   - Fit parametric model to response
   - Calculate significance based on threshold within specified time window

Returns:
    ROIAnalysis object containing processed data for each ROI
"""
function process_rois(data::Experiment{TWO_PHOTON, T}; 
    stim_indices=nothing,
    channels=nothing, 
    roi_indices=nothing,
    delay_time=50.0,
    window::Int=15,
    n_stds=2.0,
    sig_window=50.0,  # Time window in ms to look for significant responses after stimulus
    lam::T=1e4, assym::T=0.075, niter::Int=100,
    kwargs...
) where T<:Real
    
    # Get all available indices if not specified
    all_stims = getStimulusEndTime(data)
    stim_indices = isnothing(stim_indices) ? (1:length(all_stims)) : stim_indices
    channels = isnothing(channels) ? (1:size(data, 3)) : channels
    roi_list = getROImask(data) |> unique |> sort
    roi_indices = isnothing(roi_indices) ? roi_list : roi_indices

    rois = Dict{Int, Vector{ROITrace}}()

    # Create progress bar for ROIs
    p = Progress(length(roi_indices), desc="Processing ROIs: ")

    # Process each ROI
    for roi_idx in roi_indices
        traces = Vector{ROITrace}()
        
        # Process each stimulus
        for stim_idx in stim_indices
            # Calculate stimulus time windows
            stim_start_time = max(all_stims[stim_idx] - delay_time, 0.0)
            stim_end_time = if stim_idx < length(all_stims)
                min(all_stims[stim_idx+1], data.t[end])
            else
                data.t[end]
            end
            
            stim_start_index = round(Int64, (stim_start_time+1)/data.dt)
            stim_end_index = round(Int64, stim_end_time/data.dt)

            # Process each channel
            for channel_idx in channels
                # Extract ROI trace
                roi_frames = getROIarr(data, roi_idx)

                roi_trace = mean(roi_frames, dims=(1))[1, stim_start_index:stim_end_index, channel_idx]
                # Calculate dF/F and get time series
                # println(stim_start_index)
                # println(stim_end_index)
                # println(roi_trace |> size)
                pre_stim_idx = round(Int64, delay_time/data.dt)
                dFoF = baseline_trace(roi_trace; 
                    stim_frame=pre_stim_idx, 
                    window=window, 
                    lam=lam, assym=assym, niter=niter,
                )
                t_series = collect(1:length(dFoF))*data.dt
                
                # Calculate significance
                delay_idx = round(Int64, delay_time/data.dt)
                base_region = dFoF[1:delay_idx]
                pos_sig_threshold = mean(base_region) + n_stds*std(base_region)
                
                # Fit parametric model
                lb = [max(maximum(dFoF)-0.5, 0.0), 0.0, 0.0, max(delay_time-10.0, 0.0)]
                p0 = [maximum(dFoF), 100.0, 10.0, delay_time]
                ub = [maximum(dFoF)+1.0, 1000.0, 1000.0, delay_time+10.0]
                
                fit_param = fit_parametric(t_series, dFoF; lb=lb, p0=p0, ub=ub)
                
                # Check significance within specified time window
                sig_window_idx = round(Int64, sig_window/data.dt)
                sig_end_idx = min(delay_idx + sig_window_idx, length(dFoF))
                max_post_stim = maximum(dFoF[delay_idx:sig_end_idx])
                is_significant = max_post_stim > pos_sig_threshold
                
                # Create ROITrace object
                trace = ROITrace(
                    roi_trace,
                    dFoF,
                    t_series,
                    stim_start_time,
                    stim_end_time,
                    channel_idx,
                    stim_idx,  # Add stimulus index
                    fit_param,
                    is_significant
                )
                push!(traces, trace)
            end
        end
        
        rois[roi_idx] = traces
        next!(p)  # Update progress bar
    end
    
    # Create and return ROIAnalysis object
    roi_analysis = ROIAnalysis(
        rois,   
        collect(stim_indices),
        collect(channels);
        delay_time=delay_time,
        window=window,
        n_stds=n_stds,
        sig_window=sig_window,  # Add sig_window to analysis parameters
        kwargs...
    )
    
    data.HeaderDict["ROI_Analysis"] = roi_analysis # Store the analysis in the experiment's HeaderDict
    return roi_analysis
end