function process_rois(data::Experiment{TWO_PHOTON, T}; 
    
    #Parameters for the specific ROI target
    channel=1,
    roi_indices=nothing,

    #Baseline correction parameters
    window::Int=40,
    baseline_divisor_start = 1,
    baseline_divisor_end = nothing,
    linear_fill_start = nothing,
    linear_fill_end = nothing,
    
    #ROI parameters
    pos_sig_level = 2.0,
    neg_sig_level = 3.0,

    #These parameters are in terms of frames. Change for ms
    stim_frame = 50, 
    sig_threshold_std_start = 1,
    sig_threshold_std_end = 5,
    sig_threshold_mean_start = 12,
    sig_threshold_mean_end = 2,
    argmax_threshold_end = 25,
    max_dfof_end = 5,
    min_dfof_end = 100,

) where T<:Real

    roi_indices = getROImask(data) |> unique
    sig_rois = falses(length(roi_indices))
    dFoF_traces = zeros(length(roi_indices), length(data.t))

    for roi_idx in roi_indices
        println("Processing ROI $roi_idx")
        roi_arr = getROIarr(data, roi_idx)
        roi_trace = mean(roi_arr, dims = 1)[1,:, channel]
        #Adjust this a bit to get a better
        dFoF = baseline_trace(roi_trace, 
            window = window, 
            baseline_divisor_start = baseline_divisor_start, 
            baseline_divisor_end = baseline_divisor_end, 
            linear_fill_start = linear_fill_start, 
            linear_fill_end = linear_fill_end)

        #These are the parameters for the signal threshold
        sig_std = std(dFoF[sig_threshold_std_start:stim_frame-sig_threshold_std_end])
        sig_mean = mean(dFoF[stim_frame-sig_threshold_mean_start: stim_frame-sig_threshold_mean_end])
        #println("sig_mean: $sig_mean, sig_std: $sig_std")
        
        #Calculate the threshold for the signal
        sig_threshold = sig_mean + sig_std * pos_sig_level
        neg_threshold = sig_mean - sig_std * neg_sig_level
        #println("sig_threshold: $sig_threshold, neg_threshold: $neg_threshold")

        trace_argmax = argmax(dFoF[stim_frame:stim_frame+argmax_threshold_end]) + stim_frame
        
        max_dfof = mean(dFoF[trace_argmax:trace_argmax+max_dfof_end])
        min_dfof = minimum(dFoF[stim_frame:stim_frame+min_dfof_end])
        over_max = max_dfof > sig_threshold
        over_min = min_dfof > neg_threshold
        if over_max && over_min
            println("ROI $roi_idx is significant")
            sig_rois[roi_idx] = true
        end

        #println("max_dfof: $max_dfof, min_dfof: $min_dfof, over_max: $over_max, over_min: $over_min ")

        dFoF_traces[roi_idx, :] = dFoF
    
    end

    return sig_rois, dFoF_traces

end