using JSON

"""
Load parameters from JSON file and return as a dictionary.

Args:
    json_file_path: Path to the JSON file containing parameters
    
Returns:
    Dictionary containing all parameters organized by category
"""
function load_parameters(json_file_path = "parameters.json")
    # Get the directory of the current file
    current_dir = dirname(@__FILE__)
    full_path = joinpath(current_dir, json_file_path)
    
    # Read and parse the JSON file
    json_string = read(full_path, String)
    parameters = JSON.parse(json_string)
    
    return parameters
end

"""
Apply baseline correction to a trace using parameters from JSON file.
Multiple dispatch version that takes parameter filename as first argument.

Args:
    parameter_fn: Path to the JSON parameter file
    data_img: The data trace to process
    stim_frame: Frame number of the stimulus (default: 1)
    
Returns:
    Baseline corrected trace
"""
function baseline_trace(parameter_fn::String, data_img; stim_frame = 1)
    parameters = load_parameters(parameter_fn)
    baseline_params = parameters["baselining"]
    
    # Use the baseline_trace function from BaselineData.jl
    # This will be available after BaselineData.jl is loaded
    return baseline_trace(data_img,
        window = baseline_params["window"],
        stim_frame = stim_frame,
        baseline_divisor_start = baseline_params["baseline_divisor_start"],
        baseline_divisor_end = baseline_params["baseline_divisor_end"],
        linear_fill_start = baseline_params["linear_fill_start"],
        linear_fill_end = baseline_params["linear_fill_end"]
    )
end

"""
Process ROIs using parameters from JSON file.
Multiple dispatch version that takes parameter filename as first argument.

Args:
    parameter_fn: Path to the JSON parameter file
    data: Experiment data
    stim_idx: Stimulus index or indices to process (default: 1)
    kwargs: Additional arguments to pass to process_rois
    
Returns:
    Significant ROIs mask(s)
"""
function process_rois(parameter_fn::String, data; stim_idx = 1, kwargs...)
    parameters = load_parameters(parameter_fn)
    baseline_params = parameters["baselining"]
    roi_params = parameters["roi_processing"]
    
    return process_rois(data, stim_idx = stim_idx;
        # Baselining parameters
        window = baseline_params["window"],
        baseline_divisor_start = baseline_params["baseline_divisor_start"],
        baseline_divisor_end = baseline_params["baseline_divisor_end"],
        linear_fill_start = baseline_params["linear_fill_start"],
        linear_fill_end = baseline_params["linear_fill_end"],
        
        # ROI processing parameters
        pos_sig_level = roi_params["pos_sig_level"],
        neg_sig_level = roi_params["neg_sig_level"],
        sig_threshold_std_start = roi_params["sig_threshold_std_start"],
        sig_threshold_std_end = roi_params["sig_threshold_std_end"],
        sig_threshold_mean_start = roi_params["sig_threshold_mean_start"],
        sig_threshold_mean_end = roi_params["sig_threshold_mean_end"],
        argmax_threshold_end = roi_params["argmax_threshold_end"],
        max_dfof_end = roi_params["max_dfof_end"],
        min_dfof_end = roi_params["min_dfof_end"];
        kwargs...
    )
end 