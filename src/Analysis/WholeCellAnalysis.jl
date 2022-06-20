
"""
This takes the threshold of the datapoints: 
    It adds 4x the standard deviation to the mean
"""
function calculate_threshold(exp::Experiment{T}; Z::Int64 = 4) where T <: Real
    sum_data = sum(exp, dims = 2)/size(exp, 2)
	std_data = std(exp, dims = 2) * Z
    return sum_data .+ std_data
end

"""
This function returns all the time stamps in a spike or burst array
    The same function exists in RetinalChaos
"""

function get_timestamps(spike_array::BitVector, tseries::Vector{T}) where T <: Real
    diff_vals = map(i -> (spike_array[i]-spike_array[i+1]), 1:length(spike_array)-1)
    #println(diff_vals)
    diff_starts = findall(x -> x==-1, diff_vals) #This is a list of all the starting points in the array
    diff_ends = findall(x -> x==1, diff_vals) #This is a list of all the ending points in the array
    #If there is one more end than start than we have to remove the last point from the diff_starts

    if spike_array[1] #This means we start out in a spike and will most likely end o
        #println("We started out in a spike, the first value will be an end spike")
        diff_ends = diff_ends[2:end]
    end
    if length(diff_starts) > length(diff_ends)  #This happens because an end point was cutoff
        diff_starts = diff_starts[1:length(diff_ends)]
    elseif length(diff_starts) < length(diff_ends) #This happens because a start point was cutoff
        diff_ends = diff_ends[2:end]
    end

    return hcat(tseries[diff_starts], tseries[diff_ends])
end

function get_timestamps(spike_array::BitArray{N}, timestamps::Vector{T}) where N where T <: Real
    tstamps = Vector{Matrix{T}}(undef, size(spike_array,1))
    for i in 1:size(spike_array, 1)
        tstamps[i] = get_timestamps(spike_array[i, :], timestamps)
    end
    return tstamps
end

# For if no range has been provided but a threshold has
function get_timestamps(exp::Experiment{T}, threshold::AbstractArray{T}, rng::Tuple; dt = :exp) where T <: Real
    if dt == :exp
        dt = exp.dt
    end
    tseries = collect(rng[1]:dt:rng[2]-dt)
    tidxs = round.(Int64, (tseries./dt).+1)
    spike_array = exp.data_array[:, tidxs, :] .> threshold
    get_timestamps(spike_array, tseries)
end

get_timestamps(exp::Experiment{T}, threshold::AbstractArray{T}; dt = :exp) where T <: Real = get_timestamps(exp, threshold, (exp.t[1], exp.t[end]), dt = dt)
get_timestamps(exp::Experiment{T}, rng::Tuple; dt = :exp, kwargs...) where T <: Real = get_timestamps(exp, calculate_threshold(exp; kwargs...), rng, dt = dt)
get_timestamps(exp::Experiment{T}; dt = :exp, kwargs...) where T <: Real = get_timestamps(exp, calculate_threshold(exp; kwargs...), (exp.t[1], exp.t[end]), dt = dt)

function extract_interval(timestamps::Matrix{T}; 
        max_duration = 10e5, max_interval = 10e5,
        min_duration = 0.0, min_interval = 0.0
    ) where T <: Real
    durations = timestamps[:, 2] .- timestamps[:,1]
    lagged_starts = timestamps[2:end,1]
    lagged_ends = timestamps[1:end-1,2]
    intervals = lagged_starts .- lagged_ends
    return durations[min_duration .< durations .< max_duration], intervals[min_interval .< intervals .< max_interval]
end

function extract_interval(timestamp_arr::Vector{Matrix{T}}; 
        flatten = true, kwargs...
    ) where T <: Real
    if flatten
        #In this case we don't necessarily need to preserve the structure data and can collapse all entries into one
        durations = T[]
        intervals = T[]
        for idx in 1:length(timestamp_arr)
            result = extract_interval(timestamp_arr[idx]; kwargs...)
            if !isnothing(result)
                push!(durations, result[1]...)
                push!(intervals, result[2]...)
            end
        end
        return durations, intervals
    else
        println("Not implemented")
    end
end

"""
This function uses the Maximum Interval Sorting method to sort bursts in a single trace. 
    It takes in timestamps and returns the burst durations and the spikes per burst
A multiple dispatch of this function allows the max_interval to be calculated on a 3D array (x, y, and time) 
"""
function max_interval_algorithim(timestamps::Matrix{T}; 
        ISIstart::T = 0.500, ISIend::T = 0.500, IBImin::T = 1.0, DURmin::T = 0.200, SPBmin::Int64 = 4, 
        verbose = false
    ) where T <: Real
    #timestamps may be in seconds, whereas this is in milliseconds
    burst_timestamps = Tuple[]
    SPB_list = Float64[]
    if isempty(timestamps)
        if verbose >= 1
            println("No spikes detected")
        end
        return nothing
    else
        #Lets organize the spikes into intervals spikes and not spikes
        results = extract_interval(timestamps, min_duration = 0.0001)
        intervals = results[2]
        if verbose
            println("Interval durations: $(results[1])")
            println("Timestamps starts $(timestamps[:, 1])")
            println("Timestamps starts $(timestamps[:, 2])")
        end
        bursting = false
        burst_start_list = T[]
        burst_end_list = T[]
        burst_start = 0.0
        burst_end = 0.0
        SPB = 0
        idx = 1
        for i = 1:length(intervals)
            if bursting == false && intervals[i] <= ISIstart #If the cell does not start as bursting and the interval is under ISI start
                bursting = true #Begin the burst
                burst_start = timestamps[i, 1] #Record the burst start
            elseif bursting == true && intervals[i] >= ISIend || i == length(intervals) #If the cell is bursting, and the interval to the next spike is greater than ISI thresh
                bursting = false #The bursting can stop
                burst_end = timestamps[i, 2] #The burst end can be recorded
                if intervals[i] >= IBImin && (burst_end - burst_start) >= DURmin && SPB >= SPBmin #If the burst meets all the correct qualifications
                    if verbose
                        println("
                        Burst #$idx successfully added at timestamp : $burst_start -> $burst_end 
                            Duration: $(burst_end - burst_start) >  $DURmin  
                            Spikes per burst: $SPB > $SPBmin
                            IBI to burst #$(idx+1): $(intervals[i])
                            "
                            )
                    end
                    push!(burst_start_list, burst_start)
                    push!(burst_end_list, burst_end)
                    #push!(burst_timestamps, (burst_start, burst_end)) #Record it
                    push!(SPB_list, SPB)
                    SPB = 0
                    idx+=1
                elseif i == length(intervals) && (burst_end - burst_start) >= DURmin && SPB >= SPBmin
                    #a weird caveat, bursting has finished but interval has never cleared the ISIend
                    if verbose
                        println("
                        Burst #$idx successfully added at timestamp : $burst_start -> $burst_end 
                            Duration: $(burst_end - burst_start) >  $DURmin  
                            Spikes per burst: $SPB > $SPBmin
                            "
                            )
                    end
                    push!(burst_start_list, burst_start)
                    push!(burst_end_list, burst_end)
                    #push!(burst_timestamps, (burst_start, burst_end)) #Record it
                    push!(SPB_list, SPB)
                    SPB = 0
                    idx+=1
                else
                    if verbose
                        println("
                        Burst did not fit recommended qualities
                            Timestamp $idx: $burst_start -> $burst_end, 
                            DUR $idx: $(burst_end - burst_start) <  $DURmin 
                            SPB $idx: $SPB < $SPBmin
                            IBI $idx: $(intervals[i])
                            "
                        )
                    end                    
                end
            end
            if bursting == true
                SPB += 1
            end
        end
        

        if length(burst_start_list) > length(burst_end_list)
            #This algorithim usually leaves one last burst off because it has no end point. We can add this
            push!(burst_end_list, burst_start_list[end] + intervals[end])
        end        
        burst_timestamps = hcat(burst_start_list, burst_end_list)
        if isempty(burst_start_list)
            return nothing
        end
        return burst_timestamps, SPB_list
    end
end

function max_interval_algorithim(timestamp_arr::Vector{Matrix{T}}; kwargs...) where T <: Real
    #In this case we don't necessarily need to preserve the structure data and can collapse all entries into one
    bursts = Matrix{T}[]
    spd = T[]
    for idx in 1:length(timestamp_arr)
        
        result = max_interval_algorithim(timestamp_arr[idx]; kwargs...)
        if !isnothing(result)
            push!(bursts, result[1])
            push!(spd, result[2]...)
        end
    end
    return bursts, spd
end