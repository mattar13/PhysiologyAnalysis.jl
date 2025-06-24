
function process_rois(data::Experiment{TWO_PHOTON, T}; 
    
    #Parameters for the specific ROI target
    channel=2,
    roi_indices=nothing,

    #Baseline correction parameters
    window::Int=40,
    baseline_divisor_start = 20,
    baseline_divisor_end = 5,
    linear_fill_start = 5,
    linear_fill_end = 100,
    
    #ROI parameters
    pos_sig_level = 2.0,
    neg_sig_level = 3.0,

    #These parameters are in terms of frames. Change for ms
    stim_idx::Union{Int, Vector{Int}, Nothing} = 1,
    stim_frame::Union{Int, Nothing} = nothing, 
    sig_threshold_std_start = 1,
    sig_threshold_std_end = 5,
    sig_threshold_mean_start = 12,
    sig_threshold_mean_end = 2,
    argmax_threshold_end = 25,
    max_dfof_end = 25,
    min_dfof_end = 100,

) where T<:Real

    if isnothing(roi_indices)
        roi_indices = getROImask(data) |> unique
    end

    # Determine which stimulus indices to process
    if isnothing(stim_idx) && isnothing(stim_frame)
        # Process all significant stimulus indices
        stim_indices_to_process = collect(1:length(getStimulusEndTime(data)))
    elseif isnothing(stim_idx) && !isnothing(stim_frame)
        # Use single stim_frame (current behavior)
        stim_indices_to_process = [1]  # Dummy index, will use stim_frame
        use_stim_frame = true
    else
        # Use provided stim_idx (single or multiple)
        stim_indices_to_process = isa(stim_idx, Int) ? [stim_idx] : stim_idx
        use_stim_frame = false
    end

    # Initialize storage for multiple stimuli
    if length(stim_indices_to_process) > 1
        sig_rois_masks = Vector{BitVector}()
        sig_rois_indices = Vector{Vector{Int}}()
    else
        sig_rois = falses(length(roi_indices))
    end

    # Process each stimulus index
    for (i, current_stim_idx) in enumerate(stim_indices_to_process)
        if length(stim_indices_to_process) > 1
            sig_rois = falses(length(roi_indices))
        end

        # Determine stimulus frame
        if !isnothing(stim_idx) && !isnothing(stim_frame) && i == 1
            # Use provided stim_frame for first iteration
            current_stim_frame = stim_frame
        elseif !isnothing(stim_idx)
            # Calculate stim_frame from stim_idx
            t_stim = getStimulusEndTime(data)[current_stim_idx]
            current_stim_frame = round(Int, t_stim ./ data.dt)
        else
            # Use provided stim_frame
            current_stim_frame = stim_frame
        end
        println("current_stim_frame: $current_stim_frame")
        #Give warnings and make sure the indices are not out of bounds
        sig_threshold_std_start_idx = max(1, current_stim_frame - sig_threshold_std_start)
        sig_threshold_std_end_idx = min(size(data, 2), current_stim_frame - sig_threshold_std_end)
        if sig_threshold_std_start_idx == 1
            @warn "Signal Threshold Std Start exceeds the Stimulus frame"
            println("\t sig_threshold_std_start_idx: $sig_threshold_std_start, \n\t current_stim_frame: $current_stim_frame")
        end
        if sig_threshold_std_end_idx == size(data, 2)
            @warn "Signal Threshold Mean End exceeds the Stimulus frame"
            println("\t sig_threshold_std_end_idx: $sig_threshold_std_end, \n\t current_stim_frame: $current_stim_frame")
        end
        
        sig_threshold_mean_start_idx = max(1, current_stim_frame - sig_threshold_mean_start)
        sig_threshold_mean_end_idx = min(size(data, 2), current_stim_frame - sig_threshold_mean_end)
        if sig_threshold_mean_start_idx == 1
            @warn "Signal Threshold Mean Start exceeds the trace frame"
            println("\t sig_threshold_mean_start_idx: $sig_threshold_mean_start, \n\t current_stim_frame: $current_stim_frame")
        end
        if sig_threshold_mean_end_idx == size(data, 2)
            @warn "Signal Threshold Mean End exceeds the trace frame"
            println("\t sig_threshold_mean_end_idx: $sig_threshold_mean_end_idx, \n\t current_stim_frame: $current_stim_frame, \n\t trace length: $(size(data, 2))")
        end

        argmax_threshold_end_idx = min(size(data, 2), current_stim_frame + argmax_threshold_end)
        if argmax_threshold_end_idx == size(data, 2)
            @warn "Argmax Threshold End exceeds the trace frame"
            println("\t argmax_threshold_end_idx: $argmax_threshold_end_idx, \n\t current_stim_frame: $current_stim_frame, \n\t trace length: $(size(data, 2))")
        end

        max_dfof_end_idx = max(1, current_stim_frame + max_dfof_end)
        min_dfof_end_idx = min(size(data, 2), current_stim_frame + min_dfof_end)
        if max_dfof_end_idx == size(data, 2)
            @warn "Max Dfof End exceeds the trace frame"
            println("\t max_dfof_end_idx: $max_dfof_end_idx, \n\t current_stim_frame: $current_stim_frame, \n\t trace length: $(size(data, 2))")
        end
        if min_dfof_end_idx == size(data, 2)
            @warn "Min Dfof End exceeds the trace frame"
            println("\t min_dfof_end_idx: $min_dfof_end_idx, \n\t current_stim_frame: $current_stim_frame, \n\t trace length: $(size(data, 2))")
        end
        #Only give warnings once
        warning = true
        for roi_idx in roi_indices
            if roi_idx == 0
                println("\t Skipping ROI $roi_idx")
                continue
            end
            #println("Processing ROI $roi_idx")
            roi_arr = getROIarr(data, roi_idx)
            roi_trace = mean(roi_arr, dims = 1)[1,:, channel]
            #Adjust this a bit to get a better
            dFoF = baseline_trace(roi_trace, 
                stim_frame = current_stim_frame,
                window = window, 
                baseline_divisor_start = baseline_divisor_start, 
                baseline_divisor_end = baseline_divisor_end, 
                linear_fill_start = linear_fill_start, 
                linear_fill_end = linear_fill_end,
                warning = warning
            )
            #Stop giving warnings
            warning = false

            sig_std = std(dFoF[sig_threshold_std_start_idx:sig_threshold_std_end_idx])
            sig_mean = mean(dFoF[sig_threshold_mean_start_idx:sig_threshold_mean_end_idx])
            #println("sig_mean: $sig_mean, sig_std: $sig_std")
            
            #Calculate the threshold for the signal
            sig_threshold = sig_mean + sig_std * pos_sig_level
            neg_threshold = sig_mean - sig_std * neg_sig_level
            #println("sig_threshold: $sig_threshold, neg_threshold: $neg_threshold")
                        
            trace_argmax = argmax(dFoF[current_stim_frame:argmax_threshold_end_idx]) + current_stim_frame
            max_dfof = mean(dFoF[trace_argmax:max_dfof_end_idx])
            min_dfof = minimum(dFoF[current_stim_frame:min_dfof_end_idx])
            #println("max_dfof: $max_dfof, min_dfof: $min_dfof")

            over_max = max_dfof > sig_threshold
            over_min = min_dfof > neg_threshold
            if over_max && over_min
                #println("ROI $roi_idx is significant")
                sig_rois[roi_idx] = true
            end

            #println(size(dFoF))    
        end

        # Store results for this stimulus
        if length(stim_indices_to_process) > 1
            push!(sig_rois_masks, sig_rois)
            push!(sig_rois_indices, findall(sig_rois))
        end
    end

    # Store results in HeaderDict
    if length(stim_indices_to_process) > 1
        # Multiple stimuli - store as vectors
        data.HeaderDict["sig_rois_mask"] = sig_rois_masks
        data.HeaderDict["sig_rois_indices"] = sig_rois_indices
        return sig_rois_masks
    else
        # Single stimulus - store as before (but with new names)
        data.HeaderDict["sig_rois_mask"] = [sig_rois]
        data.HeaderDict["sig_rois_indices"] = [findall(sig_rois)]
        # Maintain backward compatibility with old keys
        data.HeaderDict["sig_rois_mask_segment"] = sig_rois
        data.HeaderDict["sig_rois_idxs"] = findall(sig_rois)
        return sig_rois
    end

end