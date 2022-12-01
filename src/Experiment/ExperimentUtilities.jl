"""
"""

function concat(data::Experiment{T}, data_add::Experiment{T}; mode::Symbol=:pad, position::Symbol=:post, kwargs...) where {T}
    new_data = deepcopy(data)
    if size(data, 2) > size(data_add, 2)
        #println("Original data larger $(size(data,2)) > $(size(data_add,2))")
        n_vals = abs(size(data, 2) - size(data_add, 2))
        if mode == :pad
            pad!(data_add, n_vals; position=position)
        elseif mode == :chop
            chop!(data, n_vals; position=position)
        end
    elseif size(data, 2) < size(data_add, 2)
        #println("Original data smaller $(size(data,2)) < $(size(data_add,2))")
        n_vals = abs(size(data, 2) - size(data_add, 2))
        if mode == :pad
            pad!(data, n_vals; position=position)
        elseif mode == :chop
            chop!(data_add, n_vals; position=position)
        end
    end

    push!(data, data_add)
    push!(data.stim_protocol, data_add.stim_protocol...)

    return new_data
end

function concat!(data::Experiment{T}, data_add::Experiment{T}; 
        mode::Symbol=:pad, position::Symbol=:post, verbose=false, 
        channel_mode::Symbol = :remove_extra,
        kwargs...) where {T}
    if size(data, 2) > size(data_add, 2)
        #println("Original data larger $(size(data,2)) > $(size(data_add,2))")
        n_vals = abs(size(data, 2) - size(data_add, 2))
        if mode == :pad
            pad!(data_add, n_vals; position=position)
        elseif mode == :chop
            chop!(data, n_vals; position=position)
        end
    elseif size(data, 2) < size(data_add, 2)
        #println("Original data smaller $(size(data,2)) < $(size(data_add,2))")
        n_vals = abs(size(data, 2) - size(data_add, 2))
        if mode == :pad
            pad!(data, n_vals; position=position)
        elseif mode == :chop
            chop!(data_add, n_vals; position=position)
        end
    end

    if size(data, 3) < size(data_add, 3)
        if verbose
            println("Concatenated file has too many channels")
            println(data_add.chNames)
        end
        #We want to remove channels that are hanging. 
        hanging_channels = findall((occursin.(data.chNames, data_add.chNames)).==false)
        data_dropped = drop(data_add, dim = 3, drop_idx = hanging_channels[1])
        push!(data, data_dropped)
        push!(data.stim_protocol, data_dropped.stim_protocol...)
    elseif size(data, 3) > size(data_add, 3)
        println("Original file has too many channels")
    else
        push!(data, data_add)
        push!(data.stim_protocol, data_add.stim_protocol...)
    end
end

function concat(filenames::Array{String,1}; kwargs...)
    #println("Data length is $(size(filenames, 1))")
    data = readABF(filenames[1]; average_sweeps=true, kwargs...)
    #IN this case we want to ensure that the stim_protocol is only 1 stimulus longer
    for filename in filenames[2:end]
        data_add = readABF(filename; average_sweeps=true, kwargs...)
        #println(size(data_add))
        concat!(data, data_add; kwargs...)
        #println(size(data, 1))
    end
    return data
end

concat(superfolder::String; kwargs...) = concat(parse_abf(superfolder); kwargs...)


"""
-------------------------------------------------------------------------------------
Experiment Utility

This function pads an experiment with n_add number of vals

This function has 2 versions. One of which is an inline version and modifies data
-------------------------------------------------------------------------------------
    padded_data = pad(data, n_add)
    pad!(data, n_add)

ARGS:
    data::Experiment{T} where T <: Real = The data to be padded
    n_add::Int64 = the number of datapoints added to the data array

KWARGS:
position::Symbol
    [DEFAULT, :post]
    {OPTIONS, :pre, :post}
    This specifies the chop or pad OPTIONS
        - Pre means that the chop or pad will be applied to the beginning of the data
        - Post means that the chop or pad will be applied to the end of the data

dims::Int64
    [DEFAULT, 2]
    This indicates the dimension to be padded. For the most part we will be padding
    the data which is in the 2nd dimension

val::T where T<:Real
    [DEFAULT, 0.0]
    This is the value that will be padded in the array
"""
function pad(trace::Experiment{T}, n_add::Int64; position::Symbol=:post, dims::Int64=2, val::T=0.0) where {T<:Real}
    data = deepcopy(trace)
    addon_size = collect(size(trace))
    addon_size[dims] = n_add
    addon = zeros(addon_size...)
    if position == :post
        data.data_array = [trace.data_array addon]
    elseif position == :pre
        data.data_array = [addon trace.data_array]
    end
    return data
end

function pad!(trace::Experiment{T}, n_add::Int64; position::Symbol=:post, dims::Int64=2, val::T=0.0) where {T<:Real}
    addon_size = collect(size(trace))
    addon_size[dims] = n_add
    addon = fill(val, addon_size...)
    if position == :post
        trace.data_array = [trace.data_array addon]
    elseif position == :pre
        trace.data_array = [addon trace.data_array]
    end
end

"""
-------------------------------------------------------------------------------------
Experiment Utility

This function removes n_chop values from either the front or back of the data

This function has 2 versions. One of which is an inline version and modifies data
-------------------------------------------------------------------------------------
    chopped_data = chop(data, n_chop)
    chop!(data, n_chop)

ARGS:
    data::Experiment{T} where T <: Real = The data to be chopped
    data_add::Experiment{T} where T<:Real = The data to be concatenated to
    filenames::Array{String, 1} = A array of filenames to be concatenated

KWARGS:
position::Symbol
    [DEFAULT, :post]
    {OPTIONS, :pre, :post}
    This specifies the chop or pad OPTIONS
        - Pre means that the chop or pad will be applied to the beginning of the data
        - Post means that the chop or pad will be applied to the end of the data

dims::Int64
    [DEFAULT, 2]
    This indicates the dimension to be padded. For the most part we will be padding
    the data which is in the 2nd dimension

val::T where T<:Real
    [DEFAULT, 0.0]
    This is the value that will be padded in the array
"""
function chop(trace::Experiment, n_chop::Int64; position::Symbol=:post, dims::Int64=2)
    data = copy(trace)
    resize_size = collect(size(trace))
    resize_size[dims] = (size(trace, dims) - n_chop)
    resize_size = map(x -> 1:x, resize_size)
    data.data_array = data.data_array[resize_size...]
    return data
end

function chop!(trace::Experiment, n_chop::Int64; position::Symbol=:post, dims::Int64=2)
    resize_size = collect(size(trace))
    resize_size[dims] = (size(trace, dims) - n_chop)
    resize_size = map(x -> 1:x, resize_size)
    trace.data_array = trace.data_array[resize_size...]
end

function drop!(trace::Experiment; dim=3, drop_idx=1)
    n_dims = collect(1:length(size(trace)))
    n_dims = [dim, n_dims[n_dims.!=dim]...]
    perm_data = permutedims(trace.data_array, n_dims)
    perm_data = perm_data[drop_idx.∉(1:size(trace, dim)), :, :]
    perm_data = permutedims(perm_data, sortperm(n_dims))
    trace.data_array = perm_data
end

function drop(trace::Experiment; kwargs...)
    trace_copy = copy(trace)
    drop!(trace_copy; kwargs...)
    return trace_copy
end


"""
If two experiments are being compared, then this function drops the second channel
"""
function match_channels(exp1::Experiment, exp2::Experiment)
    if size(exp1) != size(exp2)
        #we want to drop the extra channel
        match_ch = findall(exp1.chNames .== exp2.chNames)
        if size(exp1, 3) > size(exp2, 3)
            exp1 = drop(exp1, drop_idx=match_ch[1])
        else
            exp2 = drop(exp2, drop_idx=match_ch[1])
        end

    end
    return (exp1, exp2)
end

"""
This is just a convinent way to write subtraction as a function
"""
function sub_exp(exp1::Experiment, exp2::Experiment)
    if size(exp1) == size(exp2)
        data = deepcopy(exp1)
        #return a new experiment? make a new exp
        data.data_array = exp1.data_array - exp2.data_array
        return data
    else #If the channels don't match, we will automatically drop the unmatching one by default
        exp1, exp2 = match_channels(exp1, exp2)
        data = deepcopy(exp1)
        #return a new experiment? make a new exp
        data.data_array = exp1.data_array - exp2.data_array
        return data
    end
end

-(exp1::Experiment, exp2::Experiment) = sub_exp(exp1, exp2)


"""
This function is much like getindex, however this passes back the entire data instead of just the data array
"""
function getdata(trace::Experiment, sweeps, timepoints, channels::Union{String,Vector{String}})
    data = deepcopy(trace) #this copies the entire 
    data.data_array = trace[sweeps, timepoints, channels]
    data.chNames = channels
    return data
end

function getdata(trace::Experiment, sweeps, timepoints, channel::Int64; verbose=false) #I don't have an idea as to why this works differently
    if verbose
        println("$trace.here")
    end
    data = deepcopy(trace)
    data.data_array = trace[sweeps, timepoints, [channel]]
    data.chNames = [trace.chNames[channel]]
    return data
end

function getdata(trace::Experiment, sweeps, timepoints, channels::UnitRange{Int64}; verbose=false)
    if verbose
        println("here")
    end
    data = deepcopy(trace) #this copies the entire 
    data.data_array = trace[sweeps, timepoints, channels]
    data.chNames = trace.chNames[channels]
    return data
end

function getdata(trace::Experiment, sweeps, timestamps::Union{Float64,StepRangeLen{Float64}}, channels; verbose=false)
    data = deepcopy(trace) #this copies the entire 
    data.data_array = trace[sweeps, timestamps, channels]
    if length(timestamps) == 3 #This case may have happened if a full range was not provided. 
        data.t = collect(timestamps[1]:trace.dt:timestamps[end])
    else
        data.t = collect(timestamps)
    end
    return data
end

getchannel(trace::Experiment, ch_idx::Int64; verbose=false) = getdata(trace, :, :, ch_idx; verbose=verbose)

"""
This iterates through all of the channels 
"""
eachchannel(trace::Experiment; verbose=false) = Iterators.map(idx -> getchannel(trace, idx; verbose=verbose), 1:size(trace, 3))


"""
This iterates through all sweeps
"""
eachsweep(trace::Experiment) = Iterators.map(idx -> getdata(trace, idx, :, :), 1:size(trace, 1))
