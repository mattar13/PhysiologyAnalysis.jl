"""
This documentation should help you figure out what you have access to

- "dx" and "dy" are the x and y limit differentials
- "img_arr" is the raw data of all frames
- "red_zstack" is the red channel of the image array
- "grn_zstack" is the green channel of the image array
- "red_zproj" and "grn_zproj" are the red and green zprojection (taken with mean function)
- "red_img" and "grn_img" are those zprojections made into RGB{Float32} images
- "composite_img" is a composite of the red and grn images
- "red_trace" and "grn_trace" are the traces of the mean fluoresence by frame
- "dff_red_zstack" is the delta f/f (f0-f)/f taken of by a rolling mean
- "dff_grn_zstack" is the dff taken using just the total mean of the trace
- "dff_red_comp_zstack" and "dff_grn_comp_zstack" are the RGB{Float32} versions
- "dff_composite_zstack" is the composite
- "dff_red_zproj" and "dff_grn_zproj" are the zstack projections of the dff images
- "dff_red_img" and "dff_grn_img" are the image versions of the red and green z projections
- "dff_composite_img" is the composite of the red and green images
- "dff_red_trace" and "dff_grn_trace" are the trace versions
- "pks_vals" are the maxima peaks and valleys of the dff_green_trace
- "grn_row_sums" is the sum of all green section rows
- "grn_sect_arr" and "red_sect_arr" are the arrays formed as a result of the peakfinding
- "sig_rois" contains arrays of significant ROIs for each channel
- "sig_traces" is a 4D array with dimensions (rois, stims, datapoints, channels)
- "sig_tseries" contains the time points for the traces

For significant ROIs and their traces:
- data["sig_rois"][channel_idx] gives you the significant ROIs for that channel
- data["sig_traces"][roi_idx, stim_idx, :, channel_idx] gives you the trace for a specific ROI, stimulus, and channel
- data["sig_tseries"][channel_idx] gives you the time points for all traces

The sig_traces array is organized as follows:
- First dimension: ROI index
- Second dimension: Stimulus index
- Third dimension: Time points
- Fourth dimension: Channel index
"""
function open2Pdata(filename;
    #Channel parameters
    split_channel = true, 
    main_channel = :red,

    #Baselineing parameters
    grn_lam = 1e4, 
    red_lam = 1e4, 
    grn_window = 5, 
    red_window = 0,
    trunc_rng = nothing, 
    pre_event_time = 20.0, 
    post_event_time = 60.0,

    #Point to the stimulus file and adjust
    stim_filename = nothing, #New option if we want to specify a 2P stimulus dataset
    stimulus_name = "IN 2", 
    stimulus_threshold = 0.5,
    spike_train = false, 
    clock_offset = 3.0, # seconds

    #What are these again? 
    red_scale = 1.5, 
    grn_scale = 3.0, 
    section_by = :green, 
    peak_width = 40, 
    peak_min = 0.2, 
    
    #Debugging and verbosity parameters
    negative_peaks = :warn, # :nothing, :warn, or :error
    verbose = 3, # 1: basic progress, 2: timing info, 3: detailed info
)
    start_time = time()
    output = Dict{String, Any}()
    
    # Helper function for verbose logging
    function log_message(level, message, extra_info=nothing)
        if verbose >= level
            if level == 1
                println(message)
            elseif level == 2
                elapsed = time() - start_time
                println("[$elapsed s] $message")
            elseif level == 3
                elapsed = time() - start_time
                println("[$elapsed s] $message")
                if !isnothing(extra_info)
                    println("   Details: $extra_info")
                end
            end
        end
    end

    log_message(1, "Loading image...")
    experiment = readImage(filename)
    if split_channel
        deinterleave!(experiment)
        log_message(2, "Channels deinterleaved")
    end

    if !isnothing(trunc_rng)
        t_begin = trunc_rng[1]
        t_end = trunc_rng[2]
        if isnothing(trunc_rng[2])
            t_end = experiment.t[end]
        end
        if isnothing(trunc_rng[1])
            t_begin = experiment.t[1]
        end
        truncate_data!(experiment, t_begin = t_begin, t_end = t_end)
        log_message(2, "Data truncated", "Range: $t_begin to $t_end")
    end

    output["experiment"] = experiment
    output["time"] = experiment.t
    output["xlims"] = xlims = experiment.HeaderDict["xrng"]
    output["ylims"] = ylims = experiment.HeaderDict["xrng"]
    output["dt"] = experiment.dt
    output["dx"] = xlims[2]-xlims[1]
    output["dy"] = ylims[2]-ylims[1]
    log_message(1, "Image loaded")
   
    if split_channel
        log_message(2, "Processing split channels")
        
        output["img_arr"] = img_arr = get_all_frames(experiment)
        output["red_zstack"] = red_zstack = img_arr[:,:,:,2]
        output["grn_zstack"] = grn_zstack = img_arr[:,:,:,1]
        output["composite_zstack"] = RGB{Float32}.(red_zstack, grn_zstack, zeros(size(red_zstack)...))
        
        output["red_zproj"] = red_zproj = project(experiment, dims = (3))[:,:,1,2]
        output["grn_zproj"] = grn_zproj = project(experiment, dims = (3))[:,:,1,1]

        output["red_img"] = red_img = red_zproj./maximum(red_zproj)*red_scale .* RGB{Float32}(1, 0, 0)
        output["grn_img"] = grn_img = grn_zproj./maximum(grn_zproj)*grn_scale .* RGB{Float32}(0, 1, 0)
        output["composite_img"] = grn_img + red_img
        log_message(2, "Z projections completed")

        output["red_trace"] = project(experiment, dims = (1,2))[1,1,:,2]
        output["grn_trace"] = project(experiment, dims = (1,2))[1,1,:,1]
        log_message(2, "Z axis traces generated")
 
        _, output["dff_grn_trace"] = _, dff_grn_trace = baseline_trace(output["grn_trace"], window = grn_window, lam = grn_lam, niter = 100)
        _, output["dff_red_trace"] = _, dff_red_trace = baseline_trace(output["red_trace"], window = red_window, lam = red_lam, niter = 100)
        log_message(2, "Delta F/F traces calculated")
    else
        log_message(2, "Processing single channel")

        if main_channel == :red
            #╔═╡Seperate out the image array and green and red zstacks
            output["img_arr"] = output["red_zstack"] = red_zstack = img_arr = get_all_frames(experiment)
            output["grn_zstack"] = zeros(size(red_zstack))
            output["composite_zstack"] = RGB{Float32}.(red_zstack, zeros(size(red_zstack)...), zeros(size(red_zstack)...))

            output["red_zproj"] = red_zproj = project(experiment, dims = (3))[:,:,1,1]
            output["grn_zproj"] = grn_zproj = zeros(size(red_zproj))

            output["red_img"] = output["composite_img"] = red_img = red_zproj./maximum(red_zproj)*red_scale .* RGB{Float32}(1, 0, 0)
            output["grn_img"] = grn_img = zeros(size(red_img))

            output["red_trace"] = red_trace = project(experiment, dims = (1,2))[1,1,:,1]
            output["grn_trace"] = grn_trace = zeros(size(red_trace))

            _, output["dff_red_trace"] = _, dff_red_trace = baseline_trace(output["red_trace"], window = red_window, lam = red_lam, niter = 100)
            output["dff_grn_trace"] = dff_grn_trace = zeros(size(dff_red_trace))


        elseif main_channel == :grn
             #╔═╡Seperate out the image array and green and red zstacks
             output["img_arr"] = output["grn_zstack"] = grn_zstack = img_arr = get_all_frames(experiment)
             output["red_zstack"] = zeros(size(grn_zstack))
             output["composite_zstack"] = RGB{Float32}.(grn_zstack, zeros(size(grn_zstack)...), zeros(size(grn_zstack)...))
 
             output["grn_zproj"] = grn_zproj = project(experiment, dims = (3))[:,:,1,1]
             output["red_zproj"] = red_zproj = zeros(size(grn_zproj))
 
             output["grn_img"] = output["composite_img"] = grn_img = grn_zproj./maximum(grn_zproj)*grn_scale .* RGB{Float32}(1, 0, 0)
             output["red_img"] = red_img = zeros(size(grn_img))
 
             output["grn_trace"] = grn_trace = project(experiment, dims = (1,2))[1,1,:,1]
             output["red_trace"] = red_trace = zeros(size(grn_trace))
 
             _, output["dff_grn_trace"] = _, dff_grn_trace = baseline_trace(output["grn_trace"], window = grn_window, lam = grn_lam, niter = 100)
             output["dff_red_trace"] = dff_red_trace = zeros(size(dff_grn_trace))
        end
        println("Z axis traces generated")


        println("delta f/f images and traces extracted")
    end

    #%% Do and plot wavefinding events
    #If we specify a 
    if isnothing(stim_filename)
        if section_by == :green
            log_message(1, "Finding peaks in green channel")
            output["pks"], output["vals"] = pks, vals = findmaxima(dff_grn_trace, peak_width)
        elseif section_by == :red
            log_message(1, "Finding peaks in red channel")
            output["pks"], output["vals"] = pks, vals = findmaxima(dff_red_trace, peak_width)
        end
        
        # Handle negative peaks
        if any(pks .< 0)
            msg = "Found $(count(pks .< 0)) negative peaks"
            if negative_peaks == :warn
                @warn msg
            elseif negative_peaks == :error
                error(msg)
            end
        end
        
        log_message(3, "Peak finding completed", "Found $(length(pks)) peaks")
    else
        log_message(1, "Using IC stimulus for peak detection")
        addStimulus!(experiment, stim_filename, stimulus_name; flatten_episodic = true, stimulus_threshold = stimulus_threshold)
        dataIC = readABF(stim_filename, flatten_episodic = true, stimulus_name = stimulus_name, stimulus_threshold = stimulus_threshold) #Open the IC data
        start2P = experiment.HeaderDict["FileStartDateTime"]-Second(clock_offset) #The computer clocks are off by 3 seconds
        startIC = dataIC.HeaderDict["FileStartDateTime"]
        t_offset = Millisecond(startIC - start2P).value/1000 
        time_offset!(dataIC, t_offset)
        stim_protocol = getStimulusProtocol(dataIC)
        spike_train_group!(stim_protocol, 3.0) #We only need to do this if there are spike trains
        
        output["dataIC"] = dataIC
        addStimulus!(experiment, dataIC, stimulus_name)
        
        if spike_train
            spike_train_group!(stim_protocol, 3.0) #We only need to do this if there are spike trains
            
            spike_train_protocol = getStimulusProtocol(experiment)
            spike_train_group!(spike_train_protocol, 3.0)
        end
        output["tstamps"] = t_stamps = getStimulusEndTime(experiment)
        output["pks"] = pks = round.(Int64, (t_stamps./experiment.dt))
        
        #output["pks"] = pks = round.(Int64, (t_episodes./experiment.dt))
        # Handle negative peaks
        log_message(3, "Peak finding completed", "Found $(length(pks)) peaks")
        if any(pks .< 0)
            msg = "Found $(count(pks .< 0)) negative peaks"
            if negative_peaks == :warn
                @warn msg
            elseif negative_peaks == :error
                error(msg)
            end
        end
    end

    pre_event_length = floor(Int64, pre_event_time/experiment.dt)
    post_event_length = floor(Int64, post_event_time/experiment.dt)
    output["sect_time"] = new_t = LinRange(-pre_event_time, post_event_time, pre_event_length+post_event_length)
    red_sect_arr = zeros(length(new_t), length(pks))
    grn_sect_arr = zeros(length(new_t), length(pks))
    for (i, pk) in enumerate(pks)
        #println(pk-pre_event_length)
        if pk-pre_event_length < 0
            idx_start = 1
        else
            idx_start = round(Int64, pk-pre_event_length)
        end
        
        #println(pk+post_event_length)
        if pk+post_event_length > length(dff_red_trace)
            idx_end = length(dff_grn_trace) 
        else
            idx_end = round(Int64, pk+post_event_length-1)
        end
        
        idx_rng = idx_start:idx_end
        #println(idx_rng)

        grn_sect = dff_grn_trace[idx_rng]
        red_sect = dff_red_trace[idx_rng]
        red_sect_arr[1:length(red_sect), i] = red_sect
        grn_sect_arr[1:length(grn_sect), i] = grn_sect
    end
    output["grn_row_sums"] = grn_row_sums = sum(grn_sect_arr, dims = 1)[1,:]
    output["red_row_sums"] = red_row_sums = sum(red_sect_arr, dims = 1)[1,:]
    output["grn_sect_arr"] = grn_sect_arr# = grn_sect_arr[:, findall(grn_row_sums .!= 0.0)]
    output["red_sect_arr"] = red_sect_arr# = red_sect_arr[:, findall(red_row_sums .!= 0.0)]

    #We need to clean empty rows
    println("Peak finding completed")

    log_message(1, "Analysis completed")
    if verbose >= 2
        total_time = time() - start_time
        println("Total processing time: $total_time seconds")
    end
    
    return output
end

"""
    convert_to_multidim_array(nested_array)

Convert a nested array structure [channel_idx][stim_idx][roi_idx][datapoint] into a multidimensional array,
padding shorter arrays with NaN values to ensure equal length.

# Arguments
- `nested_array`: A nested array structure where the innermost arrays may have different lengths

# Returns
- A 4D array with dimensions (roi_idx, stim_idx, max_datapoints, channel_idx)
- The maximum length of datapoints found in any array
"""
function convert_to_multidim_array(nested_array)
    # Find the maximum length of datapoint arrays
    n_rois = length(nested_array[1][1])
    n_stims = length(nested_array[1])
    n_timepoints = length(nested_array[1][1][1])
    n_channels = length(nested_array)

    sig_traces = zeros(n_stims, n_rois, n_timepoints, n_channels)

    for STIM in 1:n_stims
        for CHANNEL in 1:n_channels
            for ROI in 1:n_rois
                # println("STIM: $STIM, CHANNEL: $CHANNEL, ROI: $ROI")
                trace = hcat(nested_array[CHANNEL][STIM][ROI]...)
                # println(trace |> size)
                #find an NaNs
                nan_idx = findfirst(isnan.(trace))
                # println("NaN at $nan_idx")
                sig_traces[STIM, ROI, :, CHANNEL] = trace
            end
        end
    end
    
    return sig_traces
end 

#Load some convienance functions
"""
    load_and_process_data(img_fn, stim_fn; 
        stimulus_name = "IN 3",
        split_channel = true, 
        main_channel = :grn, 
        post_event_time = 120.0,
        n_splits = 16,
        n_stds = 5.0,
        spike_train = true
    )

Core function for loading and processing 2-photon imaging data with automatic ROI analysis. This function is used
by both `load_puffing_data` and `load_electric_data` to handle the common data processing pipeline.

# Arguments
- `img_fn`: Path to the 2-photon imaging data file
- `stim_fn`: Path to the stimulus data file (typically an ABF file)

# Keyword Arguments
- `stimulus_name`: Name of the stimulus channel in the ABF file (default: "IN 3")
- `split_channel`: Whether to split the image data into separate channels (default: true)
- `main_channel`: Primary channel to analyze (:red or :grn) (default: :grn)
- `post_event_time`: Time in seconds to analyze after each stimulus (default: 120.0)
- `n_splits`: Number of splits for ROI analysis (default: 16)
- `n_stds`: Number of standard deviations for ROI detection threshold (default: 5.0)
- `spike_train`: Whether to process spike train data (default: true)

# Returns
A dictionary containing:
- Processed image data
- ROI analysis results
- Stimulus timing information
- Delta F/F traces
- Peak detection results
- Spike train information (if spike_train=true)

# Example
```julia
data = load_and_process_data(
    "path/to/image.tif",
    "path/to/stimulus.abf",
    stimulus_name = "IN 3",
    post_event_time = 60.0,
    n_splits = 32,
    spike_train = true
)
```

# Notes
- This function is typically not called directly, but rather through the specialized functions
  `load_puffing_data` or `load_electric_data`
- The function performs ROI analysis using pixel splitting and standard deviation-based thresholding
- All timing information is synchronized between the imaging and stimulus data
"""
function load_and_process_data(img_fn, stim_fn; 

    stimulus_name = "IN 3",
    stimulus_threshold = 0.5,
    spike_train = true,

    #Region of interest parameters 
    n_splits = 16,
    n_stds = 5.0,

    #Baselineing parameters
    grn_lam = 1e4, 
    red_lam = 1e4, 
    grn_window = 5, 
    red_window = 0,
    trunc_rng = nothing, 
    pre_event_time = 20.0, 
    post_event_time = 60.0,
    kwargs...
)
    println("Loading data from $(basename(img_fn))...")
    
    # Load the data
    data = open2Pdata(img_fn, 
        stim_filename = stim_fn, 
        stimulus_name = stimulus_name, 
        stimulus_threshold = stimulus_threshold, 
        spike_train = spike_train,
        grn_lam = grn_lam, 
        red_lam = red_lam, 
        grn_window = grn_window, 
        red_window = red_window,
        trunc_rng = trunc_rng,
        pre_event_time = pre_event_time,
        post_event_time = post_event_time,
        kwargs...
    )
    
    # Extract experiment and perform ROI analysis
    exp = data["experiment"]
    pixel_splits_roi!(exp, n_splits)
    roi_analysis = process_rois(exp; 
        n_stds = n_stds, 
        pre_event_time = pre_event_time, 
        post_event_time = post_event_time,
        grn_lam = grn_lam, 
        red_lam = red_lam, 
        grn_window = grn_window, 
        red_window = red_window,
    )
    
    # Add ROI analysis to the data dictionary
    data["roi_analysis"] = roi_analysis
    
    #println("Getting significant ROIs")
    
    # Collect all traces
    all_traces = []
    all_tseries = []
    all_sig_rois = []
    
    for channel_idx in axes(exp, 3)
        if main_channel == :grn
            sig_rois = all_sig_rois = get_significant_rois(roi_analysis, channel_idx = 1)
        elseif main_channel == :red 
            sig_rois = all_sig_rois = get_significant_rois(roi_analysis, channel_idx = 2)
        else
            println("Not really implemented, need to fix")
            sig_rois = all_sig_rois = get_significant_rois(roi_analysis, channel_idx = channel_idx)
        end
        #println("Processing channel $channel_idx")
        # First get significant ROIs for this channel
        #println("Found $(length(sig_rois)) significant ROIs")
        
        channel_traces = []
        for stim_idx in eachindex(data["pks"])
            #println("Processing stimulus $stim_idx")
            # Get traces only for significant ROIs
            traces = get_dfof_traces(roi_analysis, sig_rois, stim_idx = stim_idx, channel_idx = channel_idx)
            if !isempty(traces)
                push!(channel_traces, traces)
                if isempty(all_tseries)
                    push!(all_tseries, traces)  # Use traces array directly
                end
            end
        end
        push!(all_traces, channel_traces)
    end

    data["sig_traces"] = convert_to_multidim_array(all_traces)
    data["sig_tseries"] = all_tseries
    data["sig_rois"] = all_sig_rois  # Store significant ROIs for each channel
    data["mean_sig_trace"] = mean(data["sig_traces"], dims = (1, 2))[1,1,:,:]
    #println("Done loading and processing data")
    return data
end

"""
    load_puffing_data(img_fn, stim_fn; 
        stimulus_name = "IN 2", 
        split_channel = true,
        main_channel = :grn, 
        post_event_time = 120.0,
        n_splits = 16,
        n_stds = 5.0
    )

Load and process puffing data with automatic ROI analysis. This function is specifically designed for analyzing 
puffing experiments where a stimulus is applied to trigger calcium release events.

# Arguments
- `img_fn`: Path to the 2-photon imaging data file
- `stim_fn`: Path to the stimulus data file (typically an ABF file)

# Keyword Arguments
- `stimulus_name`: Name of the stimulus channel in the ABF file (default: "IN 2")
- `split_channel`: Whether to split the image data into separate channels (default: true)
- `main_channel`: Primary channel to analyze (:red or :grn) (default: :grn)
- `post_event_time`: Time in seconds to analyze after each stimulus (default: 120.0)
- `n_splits`: Number of splits for ROI analysis (default: 16)
- `n_stds`: Number of standard deviations for ROI detection threshold (default: 5.0)

# Returns
A dictionary containing:
- Processed image data
- ROI analysis results
- Stimulus timing information
- Delta F/F traces
- Peak detection results

# Example
```julia
data = load_puffing_data(
    "path/to/image.tif",
    "path/to/stimulus.abf",
    stimulus_name = "IN 2",
    post_event_time = 60.0
)
```
"""
function load_puffing_data(img_fn, stim_fn; 
    trunc_rng = nothing,
    stimulus_name = "IN 2", 
    split_channel = true,
    main_channel = :grn, 
    pre_event_time = 50.0,
    post_event_time = 60.0,
    n_splits = 16,
    n_stds = 5.0,
    red_lam = 1e4, red_window = 5,
    grn_lam = 1e4, grn_window = 5,
)
    return load_and_process_data(img_fn, stim_fn;
        trunc_rng = trunc_rng,
        stimulus_name = stimulus_name,
        split_channel = split_channel,
        main_channel = main_channel,
        pre_event_time = pre_event_time,
        post_event_time = post_event_time,
        n_splits = n_splits,
        n_stds = n_stds,
        spike_train = false,
        red_lam = red_lam, red_window = red_window,
        grn_lam = grn_lam, grn_window = grn_window,
    )
end

"""
    load_electric_data(img_fn, stim_fn; 
        stimulus_name = "IN 3", 
        split_channel = true,
        main_channel = :grn, 
        post_event_time = 120.0,
        n_splits = 16,
        n_stds = 5.0,
        red_lam = 1e4, red_window = 5,
        grn_lam = 1e4, grn_window = 5,
    )

Load and process electrical stimulation data with automatic ROI analysis. This function is specifically designed 
for analyzing experiments where electrical stimulation is used to trigger responses.

# Arguments
- `img_fn`: Path to the 2-photon imaging data file
- `stim_fn`: Path to the stimulus data file (typically an ABF file)

# Keyword Arguments
- `stimulus_name`: Name of the stimulus channel in the ABF file (default: "IN 3")
- `split_channel`: Whether to split the image data into separate channels (default: true)
- `main_channel`: Primary channel to analyze (:red or :grn) (default: :grn)
- `post_event_time`: Time in seconds to analyze after each stimulus (default: 120.0)
- `n_splits`: Number of splits for ROI analysis (default: 16)
- `n_stds`: Number of standard deviations for ROI detection threshold (default: 5.0)

# Returns
A dictionary containing:
- Processed image data
- ROI analysis results
- Stimulus timing information
- Delta F/F traces
- Peak detection results
- Spike train information (if present)

# Example
```julia
data = load_electric_data(
    "path/to/image.tif",
    "path/to/stimulus.abf",
    stimulus_name = "IN 3",
    post_event_time = 60.0,
    n_splits = 32  # More splits for finer ROI analysis
)
```
"""
function load_electric_data(img_fn, stim_fn; 
    trunc_rng = nothing,
    stimulus_name = "IN 3", 
    split_channel = true,
    main_channel = :grn, 
    pre_event_time = 50.0,
    post_event_time = 60.0,
    n_splits = 16,
    n_stds = 5.0,
    red_lam = 1e4, red_window = 5,
    grn_lam = 1e4, grn_window = 5,
)
    return load_and_process_data(img_fn, stim_fn;
        trunc_rng = trunc_rng,
        stimulus_name = stimulus_name,
        split_channel = split_channel,
        main_channel = main_channel,
        pre_event_time = pre_event_time,
        post_event_time = post_event_time,
        n_splits = n_splits,
        n_stds = n_stds,
        spike_train = true,
        red_lam = red_lam, red_window = red_window,
        grn_lam = grn_lam, grn_window = grn_window,
    )
end 

