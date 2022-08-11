#TODO: Eventually I want to add all filter functions into a single function 

"""
NEED DOCUMENTATION
"""
function filter_data(trace::Experiment{T}; kwargs...) where {T<:Real}
    data = deepcopy(trace)
    filter_data!(data; kwargs...)
    return data
end

function filter_data!(trace::Experiment{T}; freq=300.0, freqstop = 1000.0, bandwidth = 10.0,
        mode = :Lowpass, method = :Butterworth, 
        pole=8, ripple = 10.0, attenuation = 100.0
    ) where {T<:Real}

    #Determine the filter response
    if mode == :Lowpass
        responsetype = Lowpass(freq; fs=1 / trace.dt)
    elseif mode == :Highpass
        responsetype = Highpass(freq; fs=1 / trace.dt)
    elseif mode == :Bandpass
        responsetype = Bandpass(freq, freqstop, fs = 1/trace.dt)
    elseif mode == :Bandstop
        responsetype = Bandstop(freq, freqstop, fs = 1/trace.dt)
    end
    
    #Determine the method for designing the filter
    if method == :Butterworth
        designmethod = Butterworth(pole)
    elseif method == :Chebyshev1
        designmethod = Chebyshev1(pole, ripple)
    elseif method == :Chebyshev2
        designmethod = Chebyshev2(pole, ripple)
    elseif method == :Elliptic
        designmethod = Elliptic(pole, ripple, attenuation)
    end

    if mode == :Notch
        digital_filter = iirnotch(freq, bandwidth, fs=1 / trace.dt)
    else
        digital_filter = digitalfilter(responsetype, designmethod)
    end


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

function cwt_filter!(trace::Experiment{T}; wave=cDb2, β::Int64=2, 
    period_window::Tuple{Int64,Int64} = (1, 9), 
    level_window::Tuple{T, T} = (-Inf, Inf),
    return_cwt = false
) where T <: Real
    c = wavelet(wave, β=β)
    if return_cwt
        cwt = Matrix{Matrix{T}}(undef, size(trace,1), size(trace,3))
        for swp = 1:size(trace, 1), ch = 1:size(trace, 3)
            cwt[swp, ch] = ContinuousWavelets.cwt(trace[swp, :, ch], c)
        end
        return cwt
    else
        for swp = 1:size(trace, 1), ch = 1:size(trace, 3)
            y = ContinuousWavelets.cwt(trace[swp, :, ch], c)
            reconstruct = zeros(size(y))
            #if period_window[end] == Inf
            #    period_window = 
            if !any(period_window .== -1)
                reconstruct[:, period_window[1]:period_window[2]] .= y[:, period_window[1]:period_window[2]]
                trace.data_array[swp, :, ch] .= ContinuousWavelets.icwt(reconstruct, c, PenroseDelta()) |> vec
            else
                #zero all numbers outside of the level window
                outside_window = findall(level_window[1] .< y .< level_window[2])
                reconstruct[:, period_window[1]:period_window[2]] .= y[:, period_window[1]:period_window[2]]
                trace.data_array[swp, :, ch] .= ContinuousWavelets.icwt(reconstruct, c, PenroseDelta()) |> vec
            
            end
        end
    end
end

"""

"""
function dwt_filter!(trace::Experiment; wave=WT.db4, period_window::Tuple{Int64,Int64}=(1, 8))
    #In this case we have to limit the analyis to the window of dyadic time
    #This means that we can only analyze sizes if they are equal to 2^dyadic
    dyad_n = trunc(Int64, log(2, size(trace, 2)))
    #println(2^dyad_n)
    if (length(trace.t) - 2^dyad_n) > 0
        @warn "Sampling rate is non-dyadic. Will result in $(length(trace.t) - 2^dyad_n) leftover points"
    end
    if period_window[2] > dyad_n
        println("Period Window larger than dyad")
        period_window = (period_window[1], dyad_n)
    end
    for swp = 1:size(trace, 1), ch = 1:size(trace, 3)
        x = trace[swp, 1:2^dyad_n, ch]
        xt = dwt(x, wavelet(wave), dyad_n)
        reconstruct = zeros(size(xt))
        reconstruct[2^period_window[1]:2^(period_window[2])] .= xt[2^period_window[1]:2^(period_window[2])]
        trace.data_array[swp, 1:2^dyad_n, ch] .= idwt(reconstruct, wavelet(wave))
    end
end

function dwt_filter(trace::Experiment; kwargs...)
    #In this case we have to limit the analyis to the window of dyadic time
    #This means that we can only analyze sizes if they are equal to 2^dyadic
    #We can fix this by taking a 
    data = deepcopy(trace)
    dwt_filter!(data; kwargs...)
    return data
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
        band!(data, center=reference_filter * cycle, std=bandpass)
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