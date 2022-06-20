"""
The files in path array or super folder are concatenated into a single Experiment file
- There are a few modes
    - pre_pad will add zeros at the beginning of a data array that is too short
    - post_pad will add zeros at the end of a data array that is too short
    - pre_chop will remove beginning datapoints of a data array that is too long
    - post_chop will remove end datapoints of a data array that is too long
    - auto mode will will select a mode for you
        - If a majority of arrays are longer, it will pad the shorter ones
        - If a majority of arrays are shorter, it will chop the longer ones
"""

function concat(data::Experiment{T}, data_add::Experiment{T}; mode = :pad, position = :post, kwargs...) where {T}
    new_data = deepcopy(data)
    if size(data, 2) > size(data_add, 2)
        #println("Original data larger $(size(data,2)) > $(size(data_add,2))")
        n_vals = abs(size(data, 2) - size(data_add, 2))
        if mode == :pad
            pad!(data_add, n_vals; position = position)
        elseif mode == :chop
            chop!(data, n_vals; position = position)
        end
    elseif size(data, 2) < size(data_add, 2)
        #println("Original data smaller $(size(data,2)) < $(size(data_add,2))")
        n_vals = abs(size(data, 2) - size(data_add, 2))
        if mode == :pad
            pad!(data, n_vals; position = position)
        elseif mode == :chop
            chop!(data_add, n_vals; position = position)
        end
    end

    push!(data, data_add)
    push!(data.stim_protocol, data_add.stim_protocol...)

    return new_data
end

function concat!(data::Experiment{T}, data_add::Experiment{T}; mode = :pad, position = :post, verbose = false, kwargs...) where {T}
    if size(data, 2) > size(data_add, 2)
        #println("Original data larger $(size(data,2)) > $(size(data_add,2))")
        n_vals = abs(size(data, 2) - size(data_add, 2))
        if mode == :pad
            pad!(data_add, n_vals; position = position)
        elseif mode == :chop
            chop!(data, n_vals; position = position)
        end
    elseif size(data, 2) < size(data_add, 2)
        #println("Original data smaller $(size(data,2)) < $(size(data_add,2))")
        n_vals = abs(size(data, 2) - size(data_add, 2))
        if mode == :pad
            pad!(data, n_vals; position = position)
        elseif mode == :chop
            chop!(data_add, n_vals; position = position)
        end
    end

    if size(data, 3) != size(data_add, 3) 
        if verbose
            println(size(data))
            println(size(data_add))
            println(data_add.infoDict["abfPath"])
            #We need to write a catch here to concatenate files with different numbers of channels
            println("Don't concatenate these files")
        end
    else
        push!(data, data_add)
        push!(data.stim_protocol, data_add.stim_protocol...)
    end
end

function concat(path_arr::Array{String,1}; kwargs...)
    data = readABF(path_arr[1]; average_sweeps = true, kwargs...)
    #IN this case we want to ensure that the stim_protocol is only 1 stimulus longer
    for path in path_arr[2:end]
        data_add = readABF(path; average_sweeps = true, kwargs...)
        concat!(data, data_add; kwargs...)
    end
    return data
end

concat(superfolder::String; kwargs...) = concat(parse_abf(superfolder); kwargs...)

function pad(trace::Experiment{T}, n_add::Int64; position::Symbol = :post, dims::Int64 = 2, val::T = 0.0) where {T}
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

function pad!(trace::Experiment{T}, n_add::Int64; position::Symbol = :post, dims::Int64 = 2, val::T = 0.0) where {T}
    addon_size = collect(size(trace))
    addon_size[dims] = n_add
    addon = fill(val, addon_size...)
    if position == :post
        trace.data_array = [trace.data_array addon]
    elseif position == :pre
        trace.data_array = [addon trace.data_array]
    end
end

function chop(trace::Experiment, n_chop::Int64; position::Symbol = :post, dims::Int64 = 2)
    data = copy(trace)
    resize_size = collect(size(trace))
    resize_size[dims] = (size(trace, dims) - n_chop)
    resize_size = map(x -> 1:x, resize_size)
    data.data_array = data.data_array[resize_size...]
    return data
end

function chop!(trace::Experiment, n_chop::Int64; position::Symbol = :post, dims::Int64 = 2)
    resize_size = collect(size(trace))
    resize_size[dims] = (size(trace, dims) - n_chop)
    resize_size = map(x -> 1:x, resize_size)
    trace.data_array = trace.data_array[resize_size...]
end

function drop!(trace::Experiment; dim = 3, drop_idx = 1)
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
            exp1 = drop(exp1, drop_idx = match_ch[1])
        else
            exp2 = drop(exp2, drop_idx = match_ch[1])
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

function getdata(trace::Experiment, sweeps, timepoints, channel::Int64; verbose = false) #I don't have an idea as to why this works differently
    if verbose
        println("$trace.here")
    end
    data = deepcopy(trace)
    data.data_array = trace[sweeps, timepoints, [channel]]
    data.chNames = [trace.chNames[channel]]
    return data
end

function getdata(trace::Experiment, sweeps, timepoints, channels::UnitRange{Int64}; verbose = false)
    if verbose
        println("here")
    end
    data = deepcopy(trace) #this copies the entire 
    data.data_array = trace[sweeps, timepoints, channels]
    data.chNames = trace.chNames[channels]
    return data
end

function getdata(trace::Experiment, sweeps, timestamps::Union{Float64,StepRangeLen{Float64}}, channels; verbose = false)
    data = deepcopy(trace) #this copies the entire 
    data.data_array = trace[sweeps, timestamps, channels]
    if length(timestamps) == 3 #This case may have happened if a full range was not provided. 
        data.t = collect(timestamps[1]:trace.dt:timestamps[end])
    else
        data.t = collect(timestamps)
    end
    return data
end

getchannel(trace::Experiment, ch_idx::Int64; verbose = false) = getdata(trace, :, :, ch_idx; verbose = verbose)

"""
This iterates through all of the channels 
"""
eachchannel(trace::Experiment; verbose = false) = Iterators.map(idx -> getchannel(trace, idx; verbose = verbose), 1:size(trace, 3))


"""
This iterates through all sweeps
"""
eachsweep(trace::Experiment) = Iterators.map(idx -> getdata(trace, idx, :, :), 1:size(trace, 1))
