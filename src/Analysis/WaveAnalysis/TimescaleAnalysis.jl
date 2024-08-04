"""
This function returns all the time stamps in a spike or burst array
    The same function exists in RetinalChaos
"""

function get_timestamps(tseries::Vector{T}, spike_array::BitArray{1}) where {T<:Real}
    diff_vals = map(i -> (spike_array[i] - spike_array[i+1]), 1:length(spike_array)-1)
    diff_starts = findall(x -> x == -1, diff_vals) #This is a list of all the starting points in the array
    diff_ends = findall(x -> x == 1, diff_vals) #This is a list of all the ending points in the array
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

    return zip(tseries[diff_starts], tseries[diff_ends]) |> collect
end

function get_timestamps(tseries::Vector{T}, spike_array::BitArray{3}) where {T<:Real}
    n_trials, n_data, n_channels = size(spike_array)
    tstamps = Array{Vector{Tuple{T,T}}}(undef, (size(spike_array,1), size(spike_array,3)))
    for trial in 1:n_trials #focus the analysis on the last data channel
        for channel in 1:n_channels
            tstamps[trial, channel] = get_timestamps(tseries, spike_array[trial, :, channel])
        end
    end
    return tstamps
end

function get_timestamps(exp::Experiment; Z = 4.0)
    thresholds = calculate_threshold(exp, Z = Z)
    spike_array = exp.data_array .> thresholds
    get_timestamps(exp.t, spike_array)  
end

function extract_interval(timestamps::Vector{Tuple{T,T}},     
    max_duration=10e5, max_interval=10e5,
    min_duration=0.0, min_interval=0.0
) where T <: Real
    durations = map(ts -> (ts[2]-ts[1]), timestamps)
    intervals = map(i -> timestamps[i][1] - timestamps[i-1][2], 2:length(timestamps))
    return durations[min_duration.<durations.<max_duration], intervals[min_interval.<intervals.<max_interval]
end

function extract_interval(timestamps::Matrix{Vector{Tuple{T,T}}}; kwargs...) where {T<:Real}
    durations = Array{Vector{T}}(undef, size(timestamps))
    intervals = Array{Vector{T}}(undef, size(timestamps))
    for trial in axes(timestamps,1)
        for channel in axes(timestamps, 2)
            tstamps = timestamps[trial, channel]
            if !isempty(tstamps)
                duration_I, interval_I = extract_interval(tstamps; kwargs...)
                durations[trial, channel] = duration_I
                intervals[trial, channel] = interval_I
            end
        end
    end
    return durations, intervals
end

"""
This function uses the Maximum Interval Sorting method to sort bursts in a single trace. 
    It takes in timestamps and returns the burst durations and the spikes per burst
"""
function max_interval_algorithim(timestamps::Vector{Tuple{T,T}};
    ISIstart::T=0.5, ISIend::T=0.5, IBImin::T=1.0, DURmin::T=0.05, SPBmin::Int64=4,
    verbose=false
) where {T<:Real}
    burst_timestamps = Tuple[]
    SPB_list = Float64[]

    #Lets organize the spipkes into intervals spikes and not spikes
    durations, intervals = extract_interval(timestamps)
    bursting = false
    burst_start_list = T[]
    burst_end_list = T[]
    burst_start = 0.0
    burst_end = 0.0
    SPB = 0
    idx = 1
    for i in eachindex(intervals)
        if verbose
            println("Spike $i")
            println(burst_start)
            println(timestamps[i][2])
        end
        #println("Are the spikes part of the burst? $bursting")
        #println("Are the intervals longer than the minimum? $(intervals[i] <= ISIstart)")
        if bursting == false && intervals[i] <= ISIstart #If the cell does not start as bursting and the interval is under ISI start
            bursting = true #Begin the burst
            burst_start = timestamps[i][1] #Record the burst start
            if verbose
                println("\t Interval begins")
            end
        elseif bursting == true && intervals[i] >= ISIend || i == length(intervals) #If the cell is bursting, and the interval to the next spike is greater than ISI thresh
            bursting = false #The bursting can stop
            burst_end = timestamps[i][2] #The burst end can be recorded
            if verbose
                println("Interval $(intervals[i]) >= IBImin $(IBImin) $(intervals[i] >= IBImin)")
                println("Duration $((burst_end - burst_start)) >= $DURmin $((burst_end - burst_start) >= DURmin)")
                println("Spike per burst $(SPB) >= $SPBmin $(SPB >= SPBmin)")
            end
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
                push!(SPB_list, SPB)
                SPB = 0
                idx += 1
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
                idx += 1
            else
                if verbose
                    println("Burst did not fit recommended qualities")
                    println("Timestamp $idx: $burst_start -> $burst_end")
                    println("DUR $idx: $(burst_end - burst_start) $(<)  $DURmin") 
                    println("SPB $idx: $SPB < $SPBmin")
                    println("IBI $idx: $(intervals[i])")
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
    return burst_timestamps, SPB_list
end

function max_interval_algorithim(timestamp_arr::VecOrMat{Matrix{T}}; reshape_back=true, kwargs...) where {T<:Real}
    n_sizes = size(timestamp_arr)
    n_flat = *(n_sizes...)
    bursts = Vector{Matrix{T}}(undef, n_flat)
    spd = Vector{Vector{T}}(undef, n_flat)
    if length(n_sizes) > 1
        resize_arr = reshape(timestamp_arr, n_flat)
    else
        #We don't need to reshape the array back since it is only 1D
        reshape_back = false
        resize_arr = timestamp_arr
    end
    for idx in 1:length(timestamp_arr)
        if isassigned(timestamp_arr, idx)
            result = max_interval_algorithim(resize_arr[idx]; kwargs...)
            bursts[idx] = result[1]
            spd[idx] = (result[2])
        end
    end
    if reshape_back
        return reshape(bursts, n_sizes...), reshape(spd, n_sizes...)
    else
        return bursts, spd
    end
end

function timeseries_analysis(t, vm_array;
    timestamps_only=false, Z::Float64=4.0,
    max_spike_duration::Float64=50.0, max_spike_interval=100,
    max_burst_duration::Float64=10e5, max_burst_interval=10e5,
    verbose=false
) #where {T, N}
    N = length(size(vm_array))
    #println(N)
    if verbose
        print("[$(now())]: Extracting the thresholds... ")
    end
    if N == 1
        thresholds = calculate_threshold(vm_array, Z=Z)
    elseif N == 2 || N == 3
        thresholds = calculate_threshold(vm_array, Z=Z, dims=2)
    end
    spike_array = Array(vm_array .> thresholds)
    #println(spike_array |> typeof)
    if verbose
        println("Completed")
    end
    spikes = get_timestamps(t, spike_array)
    res = max_interval_algorithim(spikes)

    if isnothing(res)
        bursts = spb = nothing
    else
        bursts, spb = res
    end

    timestamps = Dict(
        "Spikes" => spikes,
        "Bursts" => bursts
    )
    if timestamps_only
        return timestamps
    else
        if verbose
            print("[$(now())]: Extracting the Data... ")
        end
        #println(spikes)
        spike_durs, isi = extract_interval(spikes, max_duration=max_spike_duration, max_interval=max_spike_interval)
        spike_dur_avg = sum(spike_durs) / length(spike_durs)
        spike_dur_sem = std(spike_durs) / sqrt(length(spike_durs))
        isi_avg = sum(isi) / length(isi)
        isi_sem = std(isi) / sqrt(length(isi))
        if !isnothing(bursts)
            burst_durs, ibi = extract_interval(bursts, max_duration=max_burst_duration, max_interval=max_burst_interval)
            burst_dur_avg = sum(burst_durs) / length(burst_durs)
            burst_dur_sem = std(burst_durs) / sqrt(length(burst_durs))
            ibi_avg = sum(ibi) / length(ibi)
            ibi_sem = std(ibi) / sqrt(length(ibi))
        else
            burst_durs = []
            ibi = []
            burst_dur_avg = NaN
            burst_dur_sem = NaN
            ibi_avg = NaN
            ibi_sem = NaN
        end

        data = Dict(
            "Time" => t,
            "DataArray" => vm_array,
            "Thresholds" => thresholds, "SpikeDurs" => spike_durs,
            "SpikeDurAvg" => spike_dur_avg,
            "SpikeDurSEM" => spike_dur_sem, "ISIs" => isi,
            "ISIAvg" => isi_avg,
            "ISISEM" => isi_sem, "BurstDurs" => burst_durs,
            "BurstDurAvg" => burst_dur_avg,
            "BurstDurSEM" => burst_dur_sem, "IBIs" => ibi,
            "IBIAvg" => ibi_avg,
            "IBISEM" => ibi_sem,
            "SpikesPerBurst" => spb
        )
        if verbose
            println("Complete")
        end
        return timestamps, data
    end
end

#Eventually I will add a @require loop here
function timeseries_analysis(t::AbstractArray{T}, vm_array::Array{T,N}, save_file::String;
    tstamps_name="timestamps", data_name="data",
    verbose=false,
    kwargs...
) where {T,N}
    timestamps, data = timeseries_analysis(t, vm_array; kwargs...)
    if verbose
        print("[$(now())]: Saving data... ")
    end
    #Uncomment to use BSON file format
    #bson("$(save_file)\\timestamps.bson", timestamps)
    #bson("$(save_file)\\data.bson", data)
    #Uncomment to use JLD2 to save the packages
    save("$(save_file)/$(tstamps_name).jld2", timestamps)
    save("$(save_file)/$(data_name).jld2", data)
    if verbose
        println("Complete")
    end
    return timestamps, data
end

timeseries_analysis(exp::Experiment; kwargs...) = timeseries_analysis(exp.t, exp.data_array; kwargs...)
timeseries_analysis(exp::Experiment, save_file::String; kwargs...) = timeseries_analysis(exp.t, exp.data_array, save_file; kwargs...)

#=================================Import this for the experiment objects=================================#
function get_timestamps(exp::Experiment{T}, threshold::AbstractArray{T}, rng::Tuple; dt=:exp) where {T<:Real}
    if dt == :exp
        dt = exp.dt
    end
    tseries = collect(rng[1]:dt:rng[2]-dt)
    tidxs = round.(Int64, (tseries ./ dt) .+ 1)
    spike_array = exp.data_array[:, tidxs, :] .> threshold
    get_timestamps(spike_array, tseries)
end

get_timestamps(exp::Experiment{T}, threshold::AbstractArray{T}; dt=:exp) where {T<:Real} = get_timestamps(exp, threshold, (exp.t[1], exp.t[end]), dt=dt)
get_timestamps(exp::Experiment{T}, rng::Tuple; dt=:exp, kwargs...) where {T<:Real} = get_timestamps(exp, calculate_threshold(exp; kwargs...), rng, dt=dt)
get_timestamps(exp::Experiment{T}; dt=:exp, kwargs...) where {T<:Real} = get_timestamps(exp, calculate_threshold(exp; kwargs...), (exp.t[1], exp.t[end]), dt=dt)
