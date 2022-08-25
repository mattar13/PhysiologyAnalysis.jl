function data_filter!(data::Experiment;
     t_pre = 1.0, t_post = 4.0, truncate_based_on = :stimulus_beginning,
     highpass = false, EI_bandpass = 100.0, lowpass = 300.0,
     dwt_periods = false, #dwt_periods = (1,9),
     cwt_periods = false #cwt_periods = (1,9)
)
     #println(t_post)
     baseline_adjust!(data, mode = :slope)
     truncate_data!(data, t_pre = t_pre, t_post = t_post, truncate_based_on = truncate_based_on)

     #We will apply several filters consecutively
     if highpass != false
          filter_data!(data, mode = :Highpass, freq = highpass) #Highpass 0.5hz
     end


     if lowpass != false
          filter_data!(data, mode = :Lowpass, freq = lowpass) #cutout all high frequency noise
     end

     if cwt_periods != false
          data = cwt_filter(data;
               period_window = cwt_periods
          )
     end

     if dwt_periods != false
          data = dwt_filter(data;
               period_window = dwt_periods
          )
     end

     data * 1000
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