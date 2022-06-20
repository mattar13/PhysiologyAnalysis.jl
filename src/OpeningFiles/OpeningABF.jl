"""
    julia> using NeuroPhys
    julia> target_path1 = "test\\to_filter.abf"
    julia> readABF(target_path1)

"""
function readABF(::Type{T}, abf_data::Union{String,Vector{UInt8}};
    sweeps=-1,
    channels::Vector{String}=["Vm_prime", "Vm_prime4"],
    average_sweeps::Bool=false,
    stimulus_name="IN 7",  #One of the best places to store digital stimuli
    stimulus_threshold::T=2.5, #This is the normal voltage rating on digital stimuli
    warn_bad_channel=false, #This will warn if a channel is improper
    flatten_episodic::Bool=false, #If the stimulation is episodic and you want it to be continuous
    time_unit=:s, #The time unit is s, change to ms
    verbose::Bool=false
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
    #return Experiment(abfInfo, dt, t, data, ch_names, ch_units, ch_telegraph, stim_protocol_by_sweep)
end

readABF(abf_path::Union{String,Vector{UInt8}}; kwargs...) = readABF(Float64, abf_path; kwargs...)

#This function utilizes concat
function readABF(abf_folder::AbstractArray{String}; average_sweeps=false, kwargs...)
    data = concat(abf_folder; kwargs...) #In the inner loop we don't want to average the sweeps
    #Save the sweep averaging for here
    if average_sweeps
        average_sweeps!(data)
    end

    return data
end
