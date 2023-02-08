"""
This function truncates the data based on the amount of time.
    In most cases we want to truncate this data by the start of the stimulus. 
    This is because the start of the stimulus should be the same response in all experiments. (0.0) 
"""
#=
function truncate_data(trace::Experiment; t_pre=1.0, t_post=4.0, t_begin = 0.0, t_end = 100.0, truncate_based_on=:stimulus_beginning)
    dt = trace.dt
    data = deepcopy(trace)
    size_of_array = 0
    if isempty(trace.stim_protocol)
        #println("No explicit stimulus has been set")
        return data
    elseif truncate_based_on == :time_range
        #Use this if there is no stimulus, but rather you want to truncate according to time
        println("truncate new use")
        start_rng = round(Int64, t_begin / dt)
        end_rng = round(Int64, t_end / dt)
        #println(start_rng)
        #println(end_rng)
        data.data_array = trace.data_array[:, start_rng:end_rng, :]
        data.t = trace.t[start_rng:end_rng] .- trace.t[start_rng]
        return data
    else
        for swp = axes(trace, 1)
            stim_protocol = trace.stim_protocol[swp]
            #We are going to iterate through each sweep and truncate it
            if truncate_based_on == :stimulus_beginning
                #This will set the beginning of the stimulus as the truncation location
                truncate_loc = stim_protocol.index_range[1]
            elseif truncate_based_on == :stimulus_end
                #This will set the beginning of the simulus as the truncation 
                truncate_loc = stim_protocol.index_range[2]
            elseif truncate_based_on == :time_range
                truncate_loc = t_pre #Set the beginning to the 
                t_pre = 0.0 #
            end
            idxs_begin = Int(t_pre / dt)
            idxs_end = Int(t_post / dt) + 1

            stim_begin_adjust = idxs_begin + (stim_protocol.index_range[1] - truncate_loc)
            stim_end_adjust = idxs_end + (stim_protocol.index_range[2] - truncate_loc)
            data.stim_protocol[swp].index_range = (stim_begin_adjust, stim_end_adjust)

            t_begin_adjust = stim_protocol.timestamps[1] - trace.t[truncate_loc+1]
            t_end_adjust = stim_protocol.timestamps[2] - trace.t[truncate_loc+1]
            data.stim_protocol[swp].timestamps = (t_begin_adjust, t_end_adjust)

            t_start = round(Int, truncate_loc - idxs_begin) #Index of truncated start point
            t_start = t_start > 0 ? t_start : 1 #If the bounds are negative indexes then reset the bounds to index 1

            t_end = round(Int, truncate_loc + idxs_end) #Index of truncated end point
            t_end = t_end < size(trace, 2) ? t_end : size(trace, 2) #If the indexes are greater than the number of datapoints then reset the indexes to n
            if size_of_array == 0
                size_of_array = t_end - t_start
                data.data_array[swp, 1:t_end-t_start+1, :] .= trace.data_array[swp, t_start:t_end, :]
            elseif size_of_array != (t_end - t_start)
                #println("Check here")
                #println(size_of_array)
                #println(t_end - t_start)
                throw(error("Inconsistant array size"))
            else
                data.data_array[swp, 1:t_end-t_start+1, :] .= trace.data_array[swp, t_start:t_end, :]
                #println("truncated array is consistant with new array")
            end
        end
        data.data_array = trace.data_array[:, 1:size_of_array, :] #remake the array with only the truncated data
        data.t = range(-t_pre, t_post, length=size_of_array)
        return data
    end
end
=#

function truncate_data!(trace::Experiment; 
    t_pre=1.0, t_post=4.0, 
    t_begin = nothing, t_end = nothing, 
    truncate_based_on=:stimulus_beginning
)
    dt = trace.dt
    size_of_array = 0
    overrun_time = 0 #This is for if t_pre is set too far before the stimulus
    if truncate_based_on == :time_range || !isnothing(t_begin) && !isnothing(t_end)
        println("running here")
        #Use this if there is no stimulus, but rather you want to truncate according to time
        start_rng = round(Int64, t_begin / dt)
        end_rng = round(Int64, t_end / dt)
        println(start_rng)
        println(end_rng)
        #println(start_rng)
        #println(end_rng)
        trace.data_array = trace.data_array[:, start_rng:end_rng, :]
        trace.t = trace.t[start_rng:end_rng] .- trace.t[start_rng]
    elseif isempty(trace.stim_protocol)
        println("No explicit stimulus has been set")
        size_of_array = round(Int64, t_post / dt)
        trace.data_array = trace.data_array[:, 1:size_of_array, :] #remake the array with only the truncated data
        trace.t = range(0.0, t_post, length=size_of_array)
    else
        for swp = axes(trace, 1)
            stim_protocol = trace.stim_protocol[swp]
            #We are going to iterate through each sweep and truncate it
            #println(trace.stim_protocol[swp].index_range)
            if truncate_based_on == :stimulus_beginning
                #This will set the beginning of the stimulus as the truncation location
                truncate_loc = stim_protocol.index_range[1]
                t_begin_adjust = 0.0
                t_end_adjust = stim_protocol.timestamps[2] - stim_protocol.timestamps[1]
            elseif truncate_based_on == :stimulus_end
                #This will set the beginning of the simulus as the truncation 
                truncate_loc = stim_protocol.index_range[2]
                t_begin_adjust = stim_protocol.timestamps[1] - stim_protocol.timestamps[2]
                t_end_adjust = 0.0
            end
            trace.stim_protocol[swp].timestamps = (t_begin_adjust, t_end_adjust)
            
            #First lets calculate how many indexes we need before the stimulus
            needed_before = round(Int, t_pre / dt)
            needed_after = round(Int, t_post / dt)
            
            #println("We need $needed_before and $needed_after indexes before and after")
            have_before = truncate_loc
            have_after = size(trace, 2) - truncate_loc
            #println("We have $have_before and $have_after indexes before and after")

            if needed_before > have_before
                #println("Not enough indexes preceed the stimulus point")
                extra_indexes = needed_before - have_before
                overrun_time = extra_indexes * dt
                #println("t_pre goes $extra_indexes indexes too far")
                idxs_begin = 1
                stim_begin_adjust = stim_protocol.index_range[1]
            else
                #println("Enough indexes preceed the stimulus point")
                idxs_begin = truncate_loc - round(Int, t_pre / dt)
                stim_begin_adjust = round(Int, t_pre / dt) 
                println(stim_begin_adjust)
            end

            if needed_after > have_after
                #println("Not enough indexes proceed the stimulus point")
                idxs_end = size(trace, 2)
            else
                #println("Enough indexes proceed the stimulus point")
                idxs_end = truncate_loc + round(Int, t_post / dt) + 1
            end

            stim_end_adjust = stim_begin_adjust + (stim_protocol.index_range[2] - stim_protocol.index_range[1])
            idxs_end = idxs_end < size(trace, 2) ? idxs_end : size(trace, 2)
            trace.stim_protocol[swp].index_range = (stim_begin_adjust, stim_end_adjust)
            #println(trace.stim_protocol[swp])

            if size_of_array == 0
                size_of_array = idxs_end - idxs_begin
            end
            trace.data_array[swp, 1:idxs_end-idxs_begin+1, :] .= trace.data_array[swp, idxs_begin:idxs_end, :]

            #println(size_of_array)
        end
        #while testing, don't change anything
        #println(size_of_array)
        trace.data_array = trace.data_array[:, 1:size_of_array, :] #remake the array with only the truncated data
        trace.t = range(-t_pre + overrun_time, t_post, length=size_of_array)
    end
end

function truncate_data(trace::Experiment; kwargs...) 
    data = deepcopy(trace)
    truncate_data!(trace)
    return data
end

"""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Data Utility

This function a splits the data by channel or sweep. 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    [chs] = split_data(data)

ARGS:
    data::Experiment{T} where T <: Real = The data to be split

KWARGS:
split_by::Symbol
    [DEFAULT, :channel]
    {[}OPTIONS, ]    
"""
function split_data(exp::Experiment; split_by::Symbol=:channel)
    if split_by == :channel
        split_exp = Array{Experiment}([])
        for ch = axes(exp, 3)
            new_data = deepcopy(exp)
            split_trace = reshape(exp[:, :, ch], (size(exp, 1), size(exp, 2), 1))
            println(size(split_trace))
            new_data.data_array = split_trace
            new_data.chNames = [exp.chNames[ch]]
            new_data.chUnits = [exp.chUnits[ch]]
            push!(split_exp, new_data)
        end
        split_exp
    elseif split_by == :sweep
        split_exp = Array{Experiment}([])
        for ch = axes(exp, 1)
            new_data = deepcopy(exp)
            split_trace = reshape(exp[:, :, ch], (size(exp, 1), size(exp, 2), 1))
            println(size(split_trace))
            new_data.data_array = split_trace
            new_data.chNames = [exp.chNames[ch]]
            new_data.chUnits = [exp.chUnits[ch]]
            push!(split_exp, new_data)
        end
        split_exp
    end
end

####################These functions are for filtering and adjusting the traces################
"""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Data Utility

This function adjusts the baseline, similar to how it is done in clampfit. 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    baseline_adjust!(data)
    data = baseline_adjust(data)

ARGS:
data::Experiment{T} where T <: Real = The data to be baselined

KWARGS:
mode::Symbol
    [DEFAULT, :slope]
    {OPTIONS, :mean, :slope}
    Specifies how the baseline should be set. 
    - Mean will subtract a single value from all points in the trace
    - Slope will fit the data to a polynomial line. 
    
polyN::Int64
    [DEFAULT, 1]
    If polynomial fit is chosen, the fit line will be a polynomial of level polyN
    1 is a linear fit (y = mx+b), 2 is a bilinear (y = ax + bx^2 + c)

reion::
    [DEFAULT]
    {OPTIONS, :prestim, :whole}
    - Prestim -> measures all time before the stimulus
    - Whole -> measures the entire trace
    - (start, end) -> a custom region

"""
function baseline_adjust(trace::Experiment{T}; kwargs...) where {T<:Real}
    data = deepcopy(trace)
    baseline_adjust!(data; kwargs...)
    return data
end

function baseline_adjust!(trace::Experiment{T};
    mode::Symbol=:slope, polyN=1, region=:prestim
) where {T<:Real}
    if isempty(trace.stim_protocol)
        #println("No stim protocol exists")
    else
        for swp in axes(trace, 1)
            if isa(region, Tuple{Float64,Float64})
                rng_begin = round(Int, region[1] / trace.dt) + 1
                if region[2] > trace.t[end]
                    rng_end = length(trace.t)
                else
                    rng_end = round(Int, region[2] / trace.dt) + 1
                end
            elseif isa(region, Tuple{Int64,Int64})
                rng_begin, rng_end = region
            elseif region == :whole
                rng_begin = 1
                rng_end = length(trace)
            elseif region == :prestim
                rng_begin = 1
                rng_end = trace.stim_protocol[swp].index_range[1] #Get the first stimulus index
            end
            for ch in axes(trace, 3)
                if mode == :mean
                    if (rng_end - rng_begin) != 0
                        baseline_adjust = sum(trace.data_array[swp, rng_begin:rng_end, ch]) / (rng_end - rng_begin)
                        #Now subtract the baseline scaling value
                        trace.data_array[swp, :, ch] .= trace.data_array[swp, :, ch] .- baseline_adjust
                    else
                        #println("no pre-stimulus range exists")
                    end
                elseif mode == :slope
                    #println(rng_begin)
                    if (rng_end - rng_begin) != 0 # && rng_begin != 1
                        pfit = PN.fit(trace.t[rng_begin:rng_end], trace[swp, rng_begin:rng_end, ch], polyN)
                        #Now offset the array by the linear range
                        trace.data_array[swp, :, ch] .= trace[swp, :, ch] - pfit.(trace.t)
                    else
                        #trace.data_array[swp, :, ch] .= trace[swp, :, ch] 
                        #println("no pre-stimulus range exists")
                    end
                end
            end
        end
    end
end

exclude(A, exclusions) = A[filter(x -> !(x ∈ exclusions), eachindex(A))]

average_sweeps!(trace::Experiment{T}) where {T<:Real} = trace.data_array = sum(trace, dims=1) / size(trace, 1)

"""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Data Utility

This function averages sweeps of the data. 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    average_sweeps!(data)
    data = average_sweeps(data)

ARGS:
    data::Experiment{T} where T <: Real = The data to be averaged
 
"""
function average_sweeps(trace::Experiment{T}) where {T<:Real}
    data = deepcopy(trace)
    average_sweeps!(data)
    return data
end

function downsample(trace::Experiment{T}, sample_rate::T) where {T<:Real}
    data = deepcopy(trace)
    downsample!(data, sample_rate)
    return data
end

"""
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Data Utility

This function reduces the points and downsamples the data. (WARNING! DATA CANNOT BE UPSAMPLED) 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    downsample!(data)
    data = downsample(data)

ARGS:
    data::Experiment{T} where {T <: Real} = The data to be averaged
    sample_rate::T where {T <: Real} = The new sample rate (must be smaller then the current sample rate)


"""
function downsample!(trace::Experiment{T}, sample_rate::T) where {T<:Real}
    old_sample_rate = 1/trace.dt
    new_dt = 1 / sample_rate
    trace.dt = new_dt #set the new dt
    trace.t = trace.t[1]:new_dt:trace.t[end] #Set the new time array
    sample_idxs = 1:round(Int64, old_sample_rate/sample_rate):size(trace, 2)
    trace.data_array = trace.data_array[:, sample_idxs, :]
end

"""
Downsamples the data to the nearest dyad for DWT 
"""
function dyadic_downsample!(trace::Experiment{T}) where {T<:Real}
    n_data = length(trace.t)
    n_dyad = 2^(trunc(log2(n_data))) |> Int64
    dyad_idxs = round.(Int64, LinRange(1, length(trace.t), n_dyad)) |> collect
    trace.t = LinRange(trace.t[1], trace.t[end], n_dyad) |> collect
    trace.dt = abs(trace.t[3] - trace.t[2])
    trace.data_array = trace.data_array[:, dyad_idxs, :]
end

function dyadic_downsample(trace::Experiment{T}) where T<:Real
    data = deepcopy(trace)
    dyadic_downsample!(data)
    return data
end