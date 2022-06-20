
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