"""
This function helps us to determine sweeps and channels in a layout for plotting
"""
function layout_helper(x::Symbol, trace_size)
     if x == :channels
          return trace_size[1]
     elseif x == :sweeps
          return trace_size[2]
     end
end
layout_helper(x::Int64, trace_size) = x

#These are all just convienance functions to help select the subplot
subplot_selector(x::Int64, trace_size) = [x]
subplot_selector(x::AbstractArray{T}, trace_size) where {T<:Real} = x
subplot_selector(x::UnitRange{T}, trace_size) where {T<:Real} = x

function subplot_selector(x::Symbol, trace_size)
     if x == :sweeps
          return 1:trace_size[1]
     elseif x == :channels
          return 1:trace_size[3]
     end
end

"""
prepares a experiment to be plotted

    -Note this only plots a single channel

"""
function plot_prep(exp::Experiment; channels=1, sweeps = :all) 
     if sweeps == :all
          return (exp.t, exp.data_array[:, :, channels]')
     else
          return (exp.t, exp.data_array[sweeps, :, channels]')
     end
end