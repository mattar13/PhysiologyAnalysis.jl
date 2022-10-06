function data_filter!(data::Experiment;
     t_pre=1.0, t_post=4.0,
     avg_swp = true,
     dwt_periods=false, #dwt_periods = (1,9),
     cwt_periods=false, #cwt_periods = (1,9)
     kwargs...
)
     #Truncate first
     baseline_adjust!(data)
     truncate_data!(data, t_pre=t_pre, t_post=t_post)
     #change from mV to uV
     data * 1000.0
     if avg_swp
          average_sweeps!(data)
     end

     # This filters the data based on the settings casette
     filter_data!(data; kwargs...) #Use the extra arguments to filter

     if cwt_periods !== false
          cwt_filter!(filtered_data;
               period_window=(WT_low_val, WT_hi_val)
          )
     end

     if dwt_periods !== false
          dwt_filter!(filtered_data;
               period_window=(WT_low_val, WT_hi_val)
          )
     end
     return data
end

function data_filter(data::Experiment; kwargs...)
     data_copy = deepcopy(data)
     data_filter!(data_copy; kwargs...)
     return data_copy
end

function data_filter(data::Tuple{Experiment{T},Experiment{T}}; kwargs...) where {T<:Real}
     data_copy = deepcopy(data[1])
     data_filter!(data_copy; kwargs...)
     return data_copy
end