function process_rois(data::Experiment{TWO_PHOTON, T}; 
    stim_indices=nothing,
    #Baseline correction parameters
    baseline_divisor_start = 1,
    baseline_divisor_end = nothing,
    linear_fill_start = nothing,
    linear_fill_end = nothing,
    #ROI parameters
    channels=nothing,
    roi_indices=nothing,
    delay_time=50.0,
    window::Int=15,
    n_stds=2.0,
    sig_window=50.0,  # Time window in ms to look for significant responses after stimulus
    analysis_window_before=100.0,  # Time window in ms before stimulus to analyze
    analysis_window_after=50.0,   # Time window in ms after stimulus to analyze
) where T<:Real

    roi_indices = getROImask(data) |> unique
    for roi_idx in roi_indices
        println("Processing ROI $roi_idx")
    end



end