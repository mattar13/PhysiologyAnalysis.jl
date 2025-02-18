#Create a new array
function map_median(input::Array{T, 2}, window_length; dim = 1, mode = :sym) where T <: Real

    output = similar(input)
    mf = MedianFilter(eltype(input), window_length)
    
    for i in axes(input, dim) # run median over each row
        # re-use mf in every iteration
        if dim == 1
            running_median!(mf, @view(output[i, :]), input[i, :], mode)
        elseif dim == 2
            running_median!(mf, @view(output[:, i]), input[:, i], mode)
        end
    end
    output
end

function mapdata_median!(exp::Experiment, window_length::Int64; channel = nothing, kwargs...)
    if isnothing(channel)
        for i in axes(exp, 3)
            mapdata_median!(exp, window_length; channel = i, kwargs...)
        end
    else
        println("channel $channel selected")
        data_channel = exp.data_array[:,:,channel]
        median_filtered = map_median(data_channel, window_length; kwargs...)
        exp.data_array[:, :, channel] .= median_filtered
    end
end

function mapdata_median(exp::Experiment, window_length::Int64; channel = nothing, kwargs...)
    exp_copy = deepcopy(exp)
    mapdata_median!(exp_copy, window_length; channel = channel, kwargs...)
    exp_copy
end