####################These functions are for filtering and adjusting the traces################
"""
This function adjusts the baseline, similar to how it is done in clampfit. 
    To change the mode of the function use the keyword argument mode
        it can cancel baseline based on: 
    - :mean -> the average voltage of a region
    - :slope -> the linear slope of a region
    To choose a region use the keyword region
    - :prestim -> measures all time before the stimulus
    - :whole -> measures the entire trace
    - (start, end) -> a custom region
It catches the baseline if the stimulus is at the beginning of the 
    """
function baseline_adjust(trace::Experiment; mode::Symbol=:slope, polyN = 1, region=:prestim)
    data = deepcopy(trace)
    if isempty(trace.stim_protocol)
        #println("No Stim protocol exists")
        return data
    else
        for swp in 1:size(trace, 1)
            if isa(region, Tuple{Float64,Float64})
                rng_begin = round(Int, region[1] / trace.dt) + 1
                if region[2] > trace.t[end]
                    rng_end = length(trace.t)
                else
                    rng_end = round(Int, region[2] / trace.dt) + 1
                end
            elseif isa(region, Tuple{Int64,Int64})
                rng_begin, rng_end = region
            elseif region == :whole
                rng_begin = 1
                rng_end = length(trace)
            elseif region == :prestim
                rng_begin = 1
                rng_end = trace.stim_protocol[swp].index_range[1] #Get the first stimulus index
            end
            for ch in 1:size(trace, 3)
                if mode == :mean
                    if (rng_end - rng_begin) != 0
                        baseline_adjust = sum(trace.data_array[swp, rng_begin:rng_end, ch]) / (rng_end - rng_begin)
                        #Now subtract the baseline scaling value
                        data.data_array[swp, :, ch] .= trace.data_array[swp, :, ch] .- baseline_adjust
                    else
                        if verbose
                            #println("no pre-stimulus range exists")
                        end
                    end
                elseif mode == :slope
                    if (rng_end - rng_begin) != 0
                        pfit = Polynomials.fit(trace.t[rng_begin:rng_end], trace[swp, rng_begin:rng_end, ch], polyN)
                        #Now offset the array by the linear range
                        data.data_array[swp, :, ch] .= trace[swp, :, ch] - pfit.(trace.t)
                    else
                        #println("no pre-stimulus range exists")
                    end
                end
            end
        end
        return data
    end
end

function baseline_adjust!(trace::Experiment; mode::Symbol=:slope, polyN = 1, region=:prestim)
    if isempty(trace.stim_protocol)
        #println("No stim protocol exists")
    else
        for swp in 1:size(trace, 1)
            if isa(region, Tuple{Float64,Float64})
                rng_begin = round(Int, region[1] / trace.dt) + 1
                if region[2] > trace.t[end]
                    rng_end = length(trace.t)
                else
                    rng_end = round(Int, region[2] / trace.dt) + 1
                end
            elseif isa(region, Tuple{Int64,Int64})
                rng_begin, rng_end = region
            elseif region == :whole
                rng_begin = 1
                rng_end = length(trace)
            elseif region == :prestim
                rng_begin = 1
                rng_end = trace.stim_protocol[swp].index_range[1] #Get the first stimulus index
            end
            for ch in 1:size(trace, 3)
                if mode == :mean
                    if (rng_end - rng_begin) != 0
                        baseline_adjust = sum(trace.data_array[swp, rng_begin:rng_end, ch]) / (rng_end - rng_begin)
                        #Now subtract the baseline scaling value
                        trace.data_array[swp, :, ch] .= trace.data_array[swp, :, ch] .- baseline_adjust
                    else
                        #println("no pre-stimulus range exists")
                    end
                elseif mode == :slope
                    #println(rng_begin)
                    if (rng_end - rng_begin) != 0 # && rng_begin != 1
                        pfit = PN.fit(trace.t[rng_begin:rng_end], trace[swp, rng_begin:rng_end, ch], polyN)
                        #Now offset the array by the linear range
                        trace.data_array[swp, :, ch] .= trace[swp, :, ch] - pfit.(trace.t)
                    else
                        #trace.data_array[swp, :, ch] .= trace[swp, :, ch] 
                        #println("no pre-stimulus range exists")
                    end
                end
            end
        end
    end
end

#TODO: Eventually I want to add all filter functions into a single function 

"""
This function applies a n-pole lowpass filter
"""
function lowpass_filter(trace::Experiment; freq=40.0, pole=8)

    responsetype = Lowpass(freq; fs=1 / trace.dt)
    designmethod = Butterworth(8)
    digital_filter = digitalfilter(responsetype, designmethod)
    data = deepcopy(trace)
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace, 3)
            #never adjust the stim
            data.data_array[swp, :, ch] .= filt(digital_filter, trace[swp, :, ch])
        end
    end
    return data
end

function lowpass_filter!(trace::Experiment; freq=40.0, pole=8)

    responsetype = Lowpass(freq; fs=1 / trace.dt)
    designmethod = Butterworth(pole)
    digital_filter = digitalfilter(responsetype, designmethod)
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace, 3)
            trace.data_array[swp, :, ch] .= filt(digital_filter, trace[swp, :, ch])
        end
    end
end

lowpass_filter(trace::Experiment, freq; pole=8) = lowpass_filter(trace; freq=freq, pole=pole)

function highpass_filter(trace::Experiment; freq=0.01, pole=8)

    responsetype = Highpass(freq; fs=1 / trace.dt)
    designmethod = Butterworth(8)
    digital_filter = digitalfilter(responsetype, designmethod)
    data = deepcopy(trace)
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace, 3)
            #never adjust the stim
            data.data_array[swp, :, ch] .= filt(digital_filter, trace[swp, :, ch])
        end
    end
    return data
end

function highpass_filter!(trace::Experiment; freq=0.01, pole=8)

    responsetype = Highpass(freq; fs=1 / trace.dt)
    designmethod = Butterworth(pole)
    digital_filter = digitalfilter(responsetype, designmethod)
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace, 3)
            trace.data_array[swp, :, ch] .= filt(digital_filter, trace[swp, :, ch])
        end
    end
end

highpass_filter(trace::Experiment, freq; pole=8) = highpass_filter(trace; freq=freq, pole=pole)

function notch_filter(trace::Experiment; center=60.0, std=30.0)
    digital_filter = iirnotch(center, std, fs=1 / trace.dt)
    data = deepcopy(trace)
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace, 3)
            #never adjust the stim
            data.data_array[swp, :, ch] .= filt(digital_filter, trace[swp, :, ch])
        end
    end
    return data
end

function notch_filter!(trace::Experiment; center=60.0, std=10.0)
    digital_filter = iirnotch(center, std, fs=1 / trace.dt)
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace, 3)
            trace.data_array[swp, :, ch] .= filt(digital_filter, trace[swp, :, ch])
        end
    end
end

function cwt_filter(trace::Experiment; wave=cDb2, β=2, dual_window=NaiveDelta(), period_window::Tuple{Int64,Int64}=(1, 9))
    data = deepcopy(trace)
    for swp = 1:size(trace, 1)
        for ch = 1:size(trace, 3)
            c = wavelet(wave, β=β)
            y = ContinuousWavelets.cwt(trace[swp, :, ch], c)
            reconstruct = zeros(size(y))
            reconstruct[:, period_window[1]:period_window[2]] .= y[:, period_window[1]:period_window[2]]
            data.data_array[swp, :, ch] .= ContinuousWavelets.icwt(reconstruct, c, dual_window) |> vec
        end
    end
    data
end

function cwt_filter!(trace::Experiment; wave=cDb2, β::Int64=2, period_window::Tuple{Int64,Int64}=(1, 9))

    for swp = 1:size(trace, 1)
        for ch = 1:size(trace, 3)
            c = wavelet(wave, β=β)
            y = ContinuousWavelets.cwt(trace[swp, :, ch], c)
            reconstruct = zeros(size(y))
            reconstruct[:, period_window[1]:period_window[2]] .= y[:, period_window[1]:period_window[2]]
            trace.data_array[swp, :, ch] .= ContinuousWavelets.icwt(reconstruct, c, PenroseDelta()) |> vec
        end
    end
end

"""

"""
function dwt_filter(trace::Experiment; wave=WT.db4, period_window::Tuple{Int64,Int64}=(1, 8), direction = :bidirectional)
    #In this case we have to limit the analyis to the window of dyadic time
    #This means that we can only analyze sizes if they are equal to 2^dyadic
    #We can fix this by taking a 
    data = deepcopy(trace)
    dyad_n = trunc(Int64, log(2, size(data, 2)))
    println(2^dyad_n)
    println(length(trace.t) - 2^dyad_n+1)
    if period_window[2] > dyad_n
        println("Period Window larger than dyad")
        period_window = (period_window[1], dyad_n)
    end
    for swp = 1:size(data, 1), ch = 1:size(data, 3)
        if direction == :forward
            x = data[swp, 1:2^dyad_n, ch]
            xt = dwt(x, wavelet(wave), dyad_n)
            reconstruct = zeros(size(xt))
            reconstruct[2^period_window[1]:2^(period_window[2])] .= xt[2^period_window[1]:2^(period_window[2])]
            data.data_array[swp, 1:2^dyad_n, ch] .= idwt(reconstruct, wavelet(wave))
        elseif direction == :reverse
            start_idx = length(trace.t) - (2^dyad_n) + 1
        
            x = data[swp, start_idx:length(data), ch]
            xt = dwt(x, wavelet(wave), dyad_n)
            reconstruct = zeros(size(xt))
            reconstruct[2^period_window[1]:2^(period_window[2])] .= xt[2^period_window[1]:2^(period_window[2])]
            data.data_array[swp, start_idx:length(data), ch] .= idwt(reconstruct, wavelet(wave))
        elseif direction == :bidirectional
            #do the reconstruction in the forward direction
            x = data[swp, 1:2^dyad_n, ch]
            xt = dwt(x, wavelet(wave), dyad_n)
            reconstruct_for = zeros(size(xt))
            reconstruct_for[2^period_window[1]:2^(period_window[2])] .= xt[2^period_window[1]:2^(period_window[2])]
        
            #Do the reconstruction in the reverse direction
            start_idx = length(trace.t) - (2^dyad_n) + 1
            x = data[swp, start_idx:length(data), ch]
            xt = dwt(x, wavelet(wave), dyad_n)
            reconstruct_rev = zeros(size(xt))
            reconstruct_rev[2^period_window[1]:2^(period_window[2])] .= xt[2^period_window[1]:2^(period_window[2])]
        
            data.data_array[swp, start_idx:length(data), ch] .= idwt(reconstruct_rev, wavelet(wave)) #We want to do reverse first 
            data.data_array[swp, 1:start_idx-1, ch] .= idwt(reconstruct_for, wavelet(wave))[1:start_idx-1] #And use the forward to fill in the first chunk
        end
    end
    data
end

"""
This is from the adaptive line interface filter in the Clampfit manual

This takes notch filters at every harmonic

#Stimulus artifacts have a very specific harmonic
250, 500, 750, 1000 ... 250n
"""
function EI_filter(trace; reference_filter=60.0, bandpass=10.0, cycles=5)
    data = deepcopy(trace)
    for cycle in 1:cycles
        notch_filter!(data, center=reference_filter * cycle, std=bandpass)
    end
    return data
end

function EI_filter!(trace; reference_filter=60.0, bandpass=10.0, cycles=5)
    for cycle in 1:cycles
        notch_filter!(trace, center=reference_filter * cycle, std=bandpass)
    end
end

function normalize(trace::Experiment; rng=(-1, 0))
    data = deepcopy(trace)
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace, 3)
            data[swp, :, ch] .= (trace[swp, :, ch] ./ minimum(trace[swp, :, ch], dims=2))
        end
    end
    return data
end

function normalize!(trace::Experiment; rng=(-1, 0))
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace, 3)
            if rng[1] < 0
                trace.data_array[swp, :, ch] .= (trace[swp, :, ch] ./ minimum(trace[swp, :, ch], dims=2))
            else
                trace.data_array[swp, :, ch] .= (trace[swp, :, ch] ./ maximum(trace[swp, :, ch], dims=2))
            end
        end
    end
end

function rolling_mean(trace::Experiment; window::Int64=10)
    data = deepcopy(trace)
    for swp in 1:size(trace, 1), ch in 1:size(trace, 3)
        for i in 1:window:size(data, 2)-window
            data.data_array[swp, i, ch] = sum(data.data_array[swp, i:i+window, ch])/window
        end
    end
    return data
end


################## Check these functions because they might be deprecated #####################################
function fft_spectrum(data::Experiment)
    #FFTW filtering
    t = data.t
    dt = t[2] - t[1]
    freqs = FFTW.fftfreq(length(t), 1.0 / dt) |> fftshift
    over_0 = findall(freqs .> 0)
    n_sweep, n_data, n_ch = size(data)
    fft_data = zeros(Complex, n_sweep, n_data, n_ch)
    for swp in 1:n_sweep
        for ch in 1:n_ch
            fft_data[swp, :, ch] = fft(data.data_array[swp, :, ch]) |> fftshift
        end
    end
    return freqs[over_0], fft_data[:, over_0, :]
end

#%% a common filter function for simplification. Remember that this is an inplace version
function filter_data!(data::Experiment;
    t_pre=1.0, t_post=4.0,
    highpass=false, EI_bandpass=100.0, lowpass=300.0,
    dwt_periods=false, #dwt_periods = (1,9),
    cwt_periods=false #cwt_periods = (1,9)
)
    #println(t_post)
    truncate_data!(data, t_pre=t_pre, t_post=t_post)
    baseline_cancel!(data, mode=:slope)

    #We will apply several filters consecutively
    if highpass != false
        highpass_filter!(data, freq=highpass) #Highpass 0.5hz
    end

    if EI_bandpass != false
        EI_filter!(data, bandpass=EI_bandpass) #adaptive line interference according to Clampfit
    end

    if lowpass != false
        lowpass_filter!(data, freq=lowpass) #cutout all high frequency noise
    end

    if cwt_periods != false
        data = cwt_filter(data;
            period_window=cwt_periods
        )
    end

    if dwt_periods != false
        data = dwt_filter(data;
            period_window=dwt_periods
        )
    end

    data * 1000
    return data
end

function filter_data(data::Experiment; kwargs...)
    data_copy = deepcopy(data)
    filter_data!(data_copy; kwargs...)
    return data_copy
end

function filter_data(data::Tuple{Experiment{T},Experiment{T}}; kwargs...) where {T<:Real}
    data_copy = deepcopy(data[1])
    filter_data!(data_copy; kwargs...)
    return data_copy
end


#This should put the data right in the middle of the cone stimuli(data, t_pre = 1.0, t_post = 1.0) #This should put the data right in the middle of the cone stimuli
cone_filter(data) = filter_data(data, t_pre=1.0, t_post=1.0)

"""
These are common and typical data filtering functions
"""
function data_filter!(data::Experiment;
    t_pre=1.0, t_post=4.0, truncate_based_on=:stimulus_beginning,
    highpass=false, EI_bandpass=100.0, lowpass=300.0,
    dwt_periods=false, #dwt_periods = (1,9),
    cwt_periods=false #cwt_periods = (1,9)
)
    #println(t_post)
    truncate_data!(data, t_pre=t_pre, t_post=t_post, truncate_based_on=truncate_based_on)
    baseline_cancel!(data, mode=:slope)

    #We will apply several filters consecutively
    if highpass != false
        highpass_filter!(data, freq=highpass) #Highpass 0.5hz
    end

    if EI_bandpass != false
        EI_filter!(data, bandpass=EI_bandpass) #adaptive line interference according to Clampfit
    end

    if lowpass != false
        lowpass_filter!(data, freq=lowpass) #cutout all high frequency noise
    end

    if cwt_periods != false
        data = cwt_filter(data;
            period_window=cwt_periods
        )
    end

    if dwt_periods != false
        data = dwt_filter(data;
            period_window=dwt_periods
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
