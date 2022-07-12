function data_filter!(data::Experiment;
     t_pre = 1.0, t_post = 4.0, truncate_based_on = :stimulus_beginning,
     highpass = false, EI_bandpass = 100.0, lowpass = 300.0,
     dwt_periods = false, #dwt_periods = (1,9),
     cwt_periods = false #cwt_periods = (1,9)
)
     #println(t_post)
     truncate_data!(data, t_pre = t_pre, t_post = t_post, truncate_based_on = truncate_based_on)
     baseline_adjust!(data, mode = :slope)

     #We will apply several filters consecutively
     if highpass != false
          highpass_filter!(data, freq = highpass) #Highpass 0.5hz
     end

     if EI_bandpass != false
          EI_filter!(data, bandpass = EI_bandpass) #adaptive line interference according to Clampfit
     end

     if lowpass != false
          lowpass_filter!(data, freq = lowpass) #cutout all high frequency noise
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