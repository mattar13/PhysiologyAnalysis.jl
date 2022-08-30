
abstract type Flash{T} end

"""
A stimulus protocol contains information about the stimulus. 
    ### T <: Real The type declaration is the type of floating point numbers
    ## Options are: 
        1) type::Symbol - This is the label or type of stimulus. TODO: This will be defined in stimulus types
        2) sweep::Int64 - The sweep number the stimulus is contained on
        3) channel::Union{String, Int64} - The channel name or number which contains stimulus data
        4) index_range::Tuple{Int64, Int64} - The indexes that thew stimulus starts at and the ends at
        5) timestamps::Tuple{T, T} - The timestamps which the stimulus starts at and ends at. 
"""
mutable struct StimulusProtocol{T}
    type::Symbol
    sweep::Int64
    channel::Union{String,Int64}
    index_range::Tuple{Int64,Int64}
    timestamps::Tuple{T,T}
end

function StimulusProtocol(type::Symbol, sweep::Int64, channel::Union{Int64,String}, index_range::Tuple{T,T}, t::Vector) where {T<:Real}
    t1 = t[index_range[1]]
    t2 = t[index_range[2]]
    StimulusProtoco(type, sweep, channel, index_range, (t1, t2))
end

#Initialize an empty stimulus protocol
StimulusProtocol() = StimulusProtocol(:None, 0, 0, (0, 0), (0.0, 0.0))

"""
this function utilizes all julia to extract ABF file data
"""

function extract_stimulus(abfInfo::Dict{String,Any}; sweep::Int64=-1, stimulus_name::String="IN 7", stimulus_threshold::Float64=2.5)
    dt = abfInfo["dataSecPerPoint"]
    stimulus_waveform = getWaveform(abfInfo, stimulus_name)
    if sweep == -1 #We want to extract info about all of the stimuli vs just one
        Stimuli = StimulusProtocol[]
        for sweep in 1:size(abfInfo["data"], 1)
            idx1 = findfirst(stimulus_waveform[sweep, :] .> stimulus_threshold)
            idx2 = findlast(stimulus_waveform[sweep, :] .> stimulus_threshold)
            if !isnothing(idx1) && !isnothing(idx2)
                push!(Stimuli, StimulusProtocol(:test, sweep, stimulus_name, (idx1, idx2), (idx1 * dt, (idx2 + 1) * dt)))
            end
        end
        return Stimuli
    else
        
        idx1 = findfirst(stimulus_waveform[sweep, :] .> stimulus_threshold)
        idx2 = findlast(stimulus_waveform[sweep, :] .> stimulus_threshold)
        if !isnothing(idx1) && !isnothing(idx2)
            return StimulusProtocol(:test, sweep, stimulus_name, (idx1, idx2), (idx1 * dt, (idx2 + 1) * dt))
        else
            return StimulusProtocol()
        end
    end
end

extract_stimulus(abf_path::String; kwargs...) = extract_stimulus(readABFInfo(abf_path); kwargs...)
