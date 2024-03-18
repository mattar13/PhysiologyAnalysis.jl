"""
    calculate_threshold(vm_arr::AbstractArray; Z = 4, dims = -1)

Finds the threshold of a trace by calculating the average and then adding the 4x the standard deviation. 
If using a differential solution, make sure dt is set, otherwise the standard deviation will be unevenly sampled
"""
function calculate_threshold(x::AbstractArray{T}; Z = 4.0, dims = -1) where {T <: Real}
    if dims == -1
        return [mean(x) + Z*std(x)]
    else
        mean = mean(x, dims = dims)
        dev = Z * std(x, dims = dims)
        return mean + dev #We want these all to come out as vectors vs matrices
    end
end

calculate_threshold(data::Experiment{F, T}; kwargs...) where {F, T<:Real} = calculate_threshold(data.data_array; kwargs...)

"""
A function that takes the threshold and returns spike array
"""

function extract_spikes(data::Experiment)

end