include("ByteMaps.jl") #These functions deal with the bytemap extractions
include("Epochs.jl") #These functions deal with the Epochs
include("WaveformExtraction.jl") #This imports the bytemaps for extracting the waveforms
include("ReadHeaders.jl")
include("ReadABFInfo.jl")

#println("ABF utilites imported")

"""
==================================================================
Reading Function

This is the baseline function for reading ABF files. 
==================================================================
    experiment = readABF(type, filename; KWARGS)
    experiment = readABF(filename; KWARGS)

ARGS:
type::Type = The type in which all data will be converted to. Defaults to Float64. 
filename::String = The filename that will be read

KWARGS:
sweeps::Union{Int64,Vector{Int64}}
    [DEFAULT, -1]
    The sweeps that will be saved. By default -1 allows all sweeps to be read. However specific sweeps can be chosen

channels::Vector{String}
    [DEFAULT, ["Vm_prime", "Vm_prime4"]] 
    The channels that will be recorded. These can be specified as a string
    By default, these are set to Vm_prime and Vm_prime4 which are voltage channels 0 and 4. 

average_sweeps::Bool 
    [DEFAULT, false]
    Specifies whether or not the function will automatically average the sweeps. This can be useful in cases where multiple files are averaged

stimulus_name::Union{String, Vector{String}}
        [DEFAULT, "IN 7"]
        This channel stores information about the stimulus. Note that this is different from the digital channel


"""
function readABF(::Type{T}, abf_data::Union{String,Vector{UInt8}};
    sweeps::Union{Int64,Vector{Int64}}=-1,
    channels::Vector{String}=["Vm_prime", "Vm_prime4"],
    average_sweeps::Bool=false,
    stimulus_name::Union{String, Vector{String}, Nothing}="IN 7",  #One of the best places to store digital stimuli
    stimulus_threshold::T=2.5, #This is the normal voltage rating on digital stimuli
    warn_bad_channel=false, #This will warn if a channel is improper
    flatten_episodic::Bool=false, #If the stimulation is episodic and you want it to be continuous
    time_unit=:s, #The time unit is s, change to ms
) where {T<:Real}
    abfInfo = readABFInfo(abf_data)
    #Pull out the requested channels
    if isa(channels, Vector{String}) #If chs is a vector of channel names extract it as such
        ch_idxs = findall(ch -> ch ∈ channels, abfInfo["adcNames"])
    elseif isa(channels, Vector{Int64}) #If chs is a vector of ints
        ch_idxs = channels
    elseif channels == -1 #if chs is -1 extract all channels
        ch_idxs = headerSection["channelList"]
    end
    #Extract info for the adc names and units
    ch_names = Vector{String}(abfInfo["adcNames"][ch_idxs])
    ch_units = Vector{String}(abfInfo["adcUnits"][ch_idxs])
    ch_telegraph = Vector{T}(abfInfo["fTelegraphAdditGain"][ch_idxs])
    #we can extract the data using getWaveform from above
    if sweeps == -1 && channels == -1
        data = abfInfo["data"]
    elseif sweeps == -1 && channels != -1
        data = getWaveform(abfInfo, ch_names; warn_bad_channel=warn_bad_channel)
    elseif sweeps != -1 && channels == -1
        data = abfInfo["data"][sweeps, :, :]
    elseif sweeps != -1 && channels != -1
        data = getWaveform(abfInfo, ch_names; warn_bad_channel=warn_bad_channel)
        data = data[sweeps, :, :]
    end
    #We need to throw an error if a dimension is empty
    if any(size(data) .== 0)
        @warn begin
            "There is in issue with the channels selected. 
            Ensure you are picking one of the following channels:
            $(abfInfo["adcNames"]) 
            or
            $(abfInfo["dacNames"])
            "
        end
        throw(DimensionMismatch)
    end
    if flatten_episodic
        n_size = size(data)
        reshape_data = permutedims(data, (3, 2, 1))
        reshape_data = reshape(reshape_data, 1, n_size[3], :)
        data = permutedims(reshape_data, (1, 3, 2))
    end


    dt = abfInfo["dataSecPerPoint"]
    t = collect(0:size(data, 2)-1) .* dt #Time is usually in seconds, but works better in ms
    if time_unit == :ms
        dt *= 1000
        t .*= 1000
    end

    stim_protocol_by_sweep = StimulusProtocol{Float64}[]
    if !isnothing(stimulus_name)
        for swp = 1:size(data, 1)
            push!(stim_protocol_by_sweep, extract_stimulus(abfInfo; sweep=swp, stimulus_name=stimulus_name, stimulus_threshold=stimulus_threshold))
        end
    end
    #This section we will rework to include getting analog and digital inputs

    if average_sweeps == true
        data = sum(data, dims=1) / size(data, 1)
        stim_protocol_by_sweep = Vector{StimulusProtocol{Float64}}([stim_protocol_by_sweep[1]])
    end
    #With our new file structure we probably need to reorganize this a bit
    return Experiment(abfInfo, dt, t, data, ch_names, ch_units, ch_telegraph, stim_protocol_by_sweep)
end

readABF(abf_path::Union{String,Vector{UInt8}}; kwargs...) = readABF(Float64, abf_path; kwargs...)

#This function utilizes concat
function OLDreadABF(filenames::AbstractArray{String}; 
    average_sweeps=true, trim_or_pad = :pad, 
    sortDate = true,
    kwargs...
)
    #old method
    data_to_cat = map(fn -> readABF(fn; kwargs...), filenames)
    if sortDate
        start_times = map(dtc -> dtc.infoDict["FileStartDateTime"], data_to_cat)
        sort_idxs = sortperm(start_times)
        data_to_cat = data_to_cat[sort_idxs]
    end
    data_sizes = map(dtc -> size(dtc, 2), data_to_cat)
    channels = map(dtc -> size(dtc, 3), data_to_cat)
    if length(unique(channels)) > 1
        #println(channels)
        #println("Concatenated file has too many channels")
        #println(data_add.chNames)
        #We want to remove channels that are hanging. 
        chNames = map(e -> e.chNames, data_to_cat)
        min_ch = minimum(channels) #This is the number of channels we need to have
        diff_idx = findall(channels .!= min_ch)
        same_idx = findfirst(channels .== min_ch)
        for di in diff_idx
            hanging_channel = findall(occursin.(chNames[same_idx], chNames[di]))
            #println(hanging_channel)
            drop!(data_to_cat[di], dim = 3, drop_idx = hanging_channel[1])
        end
    end

    if trim_or_pad == :pad
        n_add = abs.(data_sizes .- maximum(data_sizes))
    elseif trim_or_pad == :trim
        n_add = abs.(data_sizes .- minimum(data_izes))
    end
        
    for (idx, data_i) in enumerate(data_to_cat)
        if average_sweeps
            average_sweeps!(data_i)
        end
        pad!(data_i, n_add[idx])
    end
    #check to ensure all arrays are similarly sized
    return vcat(data_to_cat...)
end


function readABF(filenames::AbstractArray{String}; average_sweeps = true, 
    channel_mode::Symbol = :remove_extra,
    mode::Symbol=:pad, position::Symbol=:post, 
    verbose=false, 
    kwargs...)
    #println("Data length is $(size(filenames, 1))")
    data = readABF(filenames[1]; kwargs...)
    if average_sweeps
        average_sweeps!(data)
    end
    #IN this case we want to ensure that the stim_protocol is only 1 stimulus longer
    for filename in filenames[2:end]
        #println("Data loading takes: ")
        data_add = readABF(filename; average_sweeps=true, kwargs...)
        #println(size(data_add))
        #println("Concatenation loading takes: ")
        if average_sweeps
            average_sweeps!(data_add)
        end
        push!(data, data_add; channel_mode = channel_mode, mode = mode, position = position)
        #println(size(data, 1))
    end
    return data
end

"""
==================================================================
Parsing Function

This function finds all files with the suffix .abf 
==================================================================
    [abf_files] = parseABF(filename; KWARGS)


ARGS:
type::Type = The type in which all data will be converted to. Defaults to Float64. 
filename::String = The filename that will be read

KWARGS:
extenstion::String
    [DEFAULT, ".abf"]
    The name of the extension.

"""
function parseABF(super_folder::String; extension::String=".abf")
    file_list = String[]
    for (root, dirs, files) in walkdir(super_folder)
        for file in files
            if file[end-3:end] == extension
                path = joinpath(root, file)
                push!(file_list, path)
            end
        end
    end
    file_list
end
