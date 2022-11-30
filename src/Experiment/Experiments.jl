"""
The experiment data object. 
This contains all the data for the sweep. 
    ### DataType is defined by {T}
    ## The options are: 
        1) infoDict::Dict{String,Any} -> The information object that contains all the ABF info extracted from the binary file
        2) dt::T -> The time differential for the y axis. This is the inverse of the sampling rate collected by the digitizers
        3) t::Vector{T} -> The y axis and timestamps for each data point.
        4) data_array::Array{T,3} -> The data collected by the digitization process. The data array is sized by {sweeps, datapoints, channels}
        5) chNames::Vector{String} -> The names of each channel in string format
        6) chUnits::Vector{String} -> The units the data collected is in labeled by channel
        7) chTelegraph::Vector{T} -> The gain on each channel.
        8) stim_protocol::Vector{StimulusProtocol{T}} -> The stimulus protocols for each sweep. 
"""
mutable struct Experiment{T}
    infoDict::Union{Dict{String, Any}, Vector{Dict{String,Any}}} #This can either be a single dict or multiple
    dt::T
    t::Vector{T}
    data_array::Array{T,3}
    chNames::Vector{String}
    chUnits::Vector{String}
    chTelegraph::Vector{T}
    stim_protocol::Vector{StimulusProtocol{T}}
end

import Base: +, -, *, / #Import these basic functions to help 
function +(trace::Experiment, val::Real) 
    data = deepcopy(trace)
    data.data_array = data.data_array .+ val
    return data
end

function -(trace::Experiment, val::Real)
    data = deepcopy(trace) 
    data.data_array = data.data_array .- val
    return data
end

function *(trace::Experiment, val::Real)
    data = deepcopy(trace) 
    data.data_array = data.data_array .* val
    return data
end

function /(trace::Experiment, val::Real)
    data = deepcopy(trace) 
    data.data_array = data.data_array ./ val
    return data
end

#if the value provided is different
function /(trace::Experiment{T}, vals::Matrix{T}) where {T<:Real}
    #This function has not been worked out yet
    if size(trace, 1) == size(vals, 1) && size(trace, 3) == size(vals, 2) #Sweeps and channels of divisor match
        println("Both Sweeps and channels")
    elseif size(trace, 1) == size(vals, 1) && !(size(trace, 3) == size(vals, 2))# only channels match
        println("Only sweeps match")
    elseif !(size(trace, 1) == size(vals, 1)) && size(trace, 3) == size(vals, 2)# only channels match
        println("O")
    end
end

#This is our inplace function for scaling. Division is done by multiplying by a fraction
scaleby!(data::Experiment{T}, val::T) where T <: Real = data.data_array = data.data_array .* val 

function scaleby!(data::Experiment{T}, val::Vector{T}) where T <: Real
    #if the val is the same length of the channels then we can 
    if length(val) == size(data, 3) #Scale by the channel
        scale = reshape(val, 1,1,size(data,3))
        data.data_array = data.data_array .* scale 
    elseif length(val) == size(data,1) #Scale by sweep
        scale = reshape(val, size(data,1),1,1)
        data.data_array = data.data_array .* scale
    else
        throw(DimensionMismatch("arrays could not be broadcast to a common size; experiment dimensions: $(size(data)) vs val length: $(length(val))"))
    end
end

function scaleby(data::Experiment{T}, val) where T<:Real
    data_copy = deepcopy(data)
    scaleby!(data_copy, val)
    return data_copy
end

import Base: size, length, getindex, setindex, sum, copy, maximum, minimum,  cumsum, argmin, argmax
import Statistics.std

#Extending for Experiment
size(trace::Experiment) = size(trace.data_array)
size(trace::Experiment, dim::Int64) = size(trace.data_array, dim)

length(trace::Experiment) = size(trace, 2)

#This is the basic functionality of getindex of experiments
getindex(trace::Experiment, sweeps::Union{Int64,UnitRange{Int64}}, timepoints::Union{Int64,UnitRange{Int64}}, channels::Union{Int64,UnitRange{Int64}}) = trace.data_array[sweeps, timepoints, channels]
#getindex(trace::Experiment, sweeps::StepRangeLen{Int64}, timepoints::StepRangeLen{Int64}, channels::StepRangeLen{Int64}) = trace[sweeps, timepoints, channels]

#This function allows you to enter in a timestamp and get the data value relating to it
function getindex(trace::Experiment, sweeps, timestamps::Union{Float64,StepRangeLen{Float64}}, channels)
    @assert timestamps[end] .< trace.t[end] "Highest timestamp too high"
    @assert timestamps[1] .>= trace.t[1] "Lowest timestamp too low"

    if length(timestamps) == 3 #This case may have happened if a full range was not provided. 
        timestamps = timestamps[1]:trace.dt:timestamps[end]
    end
    offset = -round(Int64, trace.t[1] / trace.dt)
    timepoints = (timestamps ./ trace.dt) .+ 1
    timepoints = (round.(Int64, timepoints))
    timepoints .+= offset
    trace[sweeps, timepoints, channels]
end

function getindex(trace::Experiment, sweeps, timestamps, channel::String)
    ch_idx = findall(trace.chNames .== channel)
    trace[sweeps, timestamps, ch_idx]
end

function getindex(trace::Experiment, sweeps, timestamps, channels::Vector{String})
    ch_idxs = map(channel -> findall(trace.chNames .== channel)[1], channels)
    trace[sweeps, timestamps, ch_idxs]
end

#Extending get index for Experiment
getindex(trace::Experiment, I...) = trace.data_array[I...]

setindex!(trace::Experiment, v, I...) = trace.data_array[I...] = v

sum(trace::Experiment; kwargs...) = sum(trace.data_array; kwargs...)

std(trace::Experiment; kwargs...) = std(trace.data_array; kwargs...)

copy(nt::Experiment) = Experiment([getfield(nt, fn) for fn in fieldnames(nt |> typeof)]...)

minimum(trace::Experiment; kwargs...) = minimum(trace.data_array; kwargs...)

maximum(trace::Experiment; kwargs...) = maximum(trace.data_array; kwargs...)

cumsum(trace::Experiment; kwargs...) = cumsum(trace.data_array; kwargs...)

argmin(trace::Experiment; dims=2) = argmin(trace.data_array, dims=dims)

argmax(trace::Experiment; dims=2) = argmax(trace.data_array, dims=dims)

"""
Adds zeros onto the end of an array
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
Removes a portion of an array
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

import Base.push!

function push!(nt::Experiment{T}, item::AbstractArray{T}; new_name="Unnamed") where {T<:Real}
   
    #All of these options assume the new data point length matches the old one
    if size(item, 2) == size(nt, 2) && size(item, 3) == size(nt, 3)
        #item = (new_sweep, datapoints, channels)
        nt.data_array = cat(nt.data_array, item, dims=1)
    elseif size(item, 1) == size(nt, 2) && size(item, 2) == size(nt, 3)
        #item = (datapoints, channels) aka a single sweep
        item = reshape(item, 1, size(item, 1), size(item, 2))
        nt.data_array = cat(nt.data_array, item, dims=1)

    elseif size(item, 1) == size(nt, 1) && size(item, 2) == size(nt, 2)
        #item = (sweeps, datapoints, new_channels) 
        nt.data_array = cat(nt.data_array, item, dims=3)
        #Because we are adding in a new channel, add the channel name
        push!(nt.chNames, new_name)

    else
        throw(error("File size incompatible with push!"))
    end
end

function push!(nt_push_to::Experiment, nt_added::Experiment)
    #push!(nt_push_to.filename, nt_added.filename...)
    push!(nt_push_to, nt_added.data_array)
    nt_push_to.infoDict = vcat(nt_push_to.infoDict, nt_added.infoDict)
    nt_push_to.stim_protocol = vcat(nt_push_to.stim_protocol, nt_added.stim_protocol)
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

import Base: reverse, reverse!

function reverse(trace::Experiment; kwargs...)
    data = deepcopy(trace)
    data.data_array = reverse(trace.data_array; kwargs...)
    return data
end

function reverse!(trace::Experiment; kwargs...)
    trace.data_array = reverse(trace.data_array; kwargs...)
end

"""
Concatenates the data. (About to be deprecated)
"""

#=
function _concat(data::Experiment{T}, data_add::Experiment{T}; mode::Symbol=:pad, position::Symbol=:post, kwargs...) where {T}
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

function _concat!(data::Experiment{T}, data_add::Experiment{T}; 
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

function _concat(filenames::Array{String,1}; kwargs...)
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

_concat(superfolder::String; kwargs...) = concat(parse_abf(superfolder); kwargs...)
=#
import Base: _cat, cat, vcat, hcat

function _cat(dims, exps::Experiment...)
    #need to check to make sure all channel sizes are equal
    new_exp = deepcopy(exps[1])
    data_arrs = map(e -> e.data_array, exps)
    new_exp.infoDict = [map(e -> e.infoDict, exps)...]
    new_exp.data_array = cat(data_arrs..., dims = dims) 
    new_exp.stim_protocol = vcat(map(e -> e.stim_protocol, exps)...)
    #println(stim_protocol |> typeof)
    if dims == 3 #alter the channels
        new_exp.chNames = vcat(map(e -> e.chNames, exps)...)
        new_exp.chUnits = vcat(map(e -> e.chUnits, exps)...)
        new_exp.chTelegraph = vcat(map(e -> e.chTelegraph, exps)...)
    end
    return new_exp  
end

cat(exps::Experiment...; dims) = _cat(dims, exps...)
vcat(exps::Experiment...) = cat(exps..., dims = 1)
hcat(exps::Experiment...) = cat(exps..., dims = 3)

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
