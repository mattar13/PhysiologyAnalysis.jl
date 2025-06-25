function findROIcentroids(array::Array{Int,2}, xlim::AbstractRange, ylim::AbstractRange)
    # Ensure the array is square
    size_x, size_y = size(array)
    @assert size_x == size_y "The input array must be square."

    # Calculate the number of elements along each axis
    n = size_x

    # Calculate the physical spacing between pixels
    dx = (xlim[2] - xlim[1])
    dy = (ylim[2] - ylim[1])

    # Find unique labels excluding 0 (assumed to be background)
    labels = setdiff(unique(array), [0])

    # Dictionary to store centroids: label => (x_centroid, y_centroid)
    centroids = Dict{Int, Tuple{Float64, Float64}}()

    for label in labels
        # Find indices where array equals the current label
        inds = findall(array .== label)

        if isempty(inds)
            continue  # Skip if no indices found for the label
        end

        # Initialize sums for x and y positions
        sum_x = 0.0
        sum_y = 0.0

        for I in inds
            i, j = Tuple(I)

            # Map array indices to physical coordinates
            x_pos = xlim[1] + (j - 0.5) * dx
            y_pos = ylim[1] + (i - 0.5) * dy

            sum_x += x_pos
            sum_y += y_pos
        end

        # Compute the centroid by averaging positions
        num_points = length(inds)
        centroid_x = sum_x / num_points
        centroid_y = sum_y / num_points

        # Store the centroid in the dictionary
        centroids[label] = (centroid_x, centroid_y)
    end

    return centroids
end

findROIcentroids(exp) = findROIcentroids(getROImask(exp), exp.HeaderDict["xrng"], exp.HeaderDict["yrng"])

"""
    segmentTraceByStimuli(t::AbstractVector, y::AbstractVector, stimulus_times::AbstractVector; 
                          time_before::Real=0.0, time_after::Real=0.0) -> Vector{Tuple{Vector, Vector}}

Segments a timeseries trace based on stimulus times with configurable time windows.

# Arguments
- `t::AbstractVector`: Time vector
- `y::AbstractVector`: Signal vector (must have same length as t)
- `stimulus_times::AbstractVector`: Vector of stimulus times
- `time_before::Real=0.0`: Time window before each stimulus (in same units as t)
- `time_after::Real=0.0`: Time window after each stimulus (in same units as t)

# Returns
- `Vector{Tuple{Vector, Vector}}`: Array of segments, where each segment is a tuple of (t_segment, y_segment)

# Example
```julia
# Segment a trace with 1 second before and 2 seconds after each stimulus
segments = segmentTraceByStimuli(time_vector, signal_vector, stimulus_times, 
                                time_before=1.0, time_after=2.0)

# Access individual segments
for (i, (t_seg, y_seg)) in enumerate(segments)
    println("Segment: points")
end
```

# Notes
- If a stimulus time is too close to the beginning or end of the trace, the segment will be truncated
- Returns empty segments for stimuli that fall outside the trace bounds
"""
function segmentTraceByStimuli(t::AbstractVector, y::AbstractVector, stimulus_times::AbstractVector; 
                              time_before::Real=20.0, time_after::Real=100.0)
    
    # Input validation
    @assert length(t) == length(y) "Time and signal vectors must have the same length"
    @assert time_before >= 0.0 "time_before must be non-negative"
    @assert time_after >= 0.0 "time_after must be non-negative"
    
    # Ensure t is sorted (required for binary search)
    if !issorted(t)
        error("Time vector must be sorted in ascending order")
    end
    
    # Calculate time step
    dt = t[2] - t[1]
    
    # Calculate the total window size in samples
    samples_before = round(Int, time_before / dt)
    samples_after = round(Int, time_after / dt)
    total_window_size = samples_before + samples_after
    
    # Initialize the aligned array
    aligned_array = zeros(eltype(y), length(stimulus_times), total_window_size)
    
    for (stim_idx, stim_time) in enumerate(stimulus_times)
        # Find the index of the stimulus time in the time vector
        stim_idx_in_trace = searchsortedfirst(t, stim_time)
        
        # If stimulus time is not found or outside bounds, skip
        if stim_idx_in_trace > length(t) || stim_idx_in_trace < 1
            continue
        end
        
        # Calculate start and end indices for the segment
        start_idx = stim_idx_in_trace - samples_before
        end_idx = stim_idx_in_trace + samples_after - 1  # -1 because we want inclusive
        
        # Handle edge cases where segment goes outside trace bounds
        if start_idx < 1
            # Segment starts before trace begins
            actual_start = 1
            array_start = 1 + (1 - start_idx)  # Offset in the aligned array
            start_idx = 1
        else
            actual_start = start_idx
            array_start = 1
        end
        
        if end_idx > length(t)
            # Segment ends after trace ends
            actual_end = length(t)
            array_end = total_window_size - (end_idx - length(t))  # Adjust end in aligned array
            end_idx = length(t)
        else
            actual_end = end_idx
            array_end = total_window_size
        end
        
        # Extract the segment from the trace
        if actual_start <= actual_end
            segment = y[actual_start:actual_end]
            
            # Place the segment in the aligned array
            aligned_array[stim_idx, array_start:array_start + length(segment) - 1] = segment
        end
    end
    
    return aligned_array
end

"""
    segmentTraceByStimuli(exp::Experiment{TWO_PHOTON, T}, stimulus_times::AbstractVector; 
                          time_before::Real=0.0, time_after::Real=0.0) where T<:Real -> Vector{Tuple{Vector, Vector}}

Segments a timeseries trace from an experiment object based on stimulus times.

# Arguments
- `exp::Experiment{TWO_PHOTON, T}`: The experiment object
- `stimulus_times::AbstractVector`: Vector of stimulus times
- `time_before::Real=0.0`: Time window before each stimulus
- `time_after::Real=0.0`: Time window after each stimulus

# Returns
- `Vector{Tuple{Vector, Vector}}`: Array of segments, where each segment is a tuple of (t_segment, y_segment)

# Example
```julia
# Get stimulus times from experiment
stimulus_times = getStimulusStartTime(experiment)

# Segment the trace with 0.5 seconds before and 1.5 seconds after each stimulus
segments = segmentTraceByStimuli(experiment, stimulus_times, 
                                time_before=0.5, time_after=1.5)
```
"""
function segmentTraceByStimuli(exp::Experiment{TWO_PHOTON, T}, stimulus_times::AbstractVector; 
                              time_before::Real=0.0, time_after::Real=0.0) where T<:Real
    
    # Use the first channel for segmentation (you can modify this as needed)
    y = exp.data_array[:, 1, 1]  # First channel, first trial
    
    return segmentTraceByStimuli(exp.t, y, stimulus_times, 
                                time_before=time_before, time_after=time_after)
end