function data_filter!(data::Experiment;
     t_pre=1.0, t_post=4.0,
     avg_swp=false,
     scale=1000.0,
     remove_global_drift = :polyfit,
     dwt_periods=false, #dwt_periods = (1,9),
     cwt_periods=false, #cwt_periods = (1,9)
<<<<<<< Updated upstream
=======
     channels = -1,
>>>>>>> Stashed changes
     kwargs...
)
     #Truncate first
     truncate_data!(data, t_pre=t_pre, t_post=t_post)
     baseline_adjust!(data)

     #Add an extra baseline step to remove global drift
     if remove_global_drift == :polyfit
          baseline_adjust!(data, polyN = 2, region = :whole)
     elseif remove_global_drift == :Highpass
          filter_data!(data, mode = :Highpass, freq_start = 0.05)
     end

     #change from mV to uV
     scaleby!(data, scale) #scale the data by the scale number (usually is conversion from mV to μV
     if avg_swp
          average_sweeps!(data)
     end
     # This filters the data based on the settings casette
     #println(maximum(data))
     #println(minimum(data))
     filter_data!(data; kwargs...) #Use the extra arguments to filter
     #println(maximum(data))
     #println(minimum(data))

<<<<<<< Updated upstream
     if cwt_periods !== false
          cwt_filter!(filtered_data;
               period_window=(WT_low_val, WT_hi_val)
          )
     end

     if dwt_periods !== false
          dwt_filter!(filtered_data;
               period_window=(WT_low_val, WT_hi_val)
          )
=======
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
     else
          #turn the channels into a list of channels
>>>>>>> Stashed changes
     end
     #return data
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