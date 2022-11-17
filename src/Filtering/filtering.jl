#TODO: Eventually I want to add all filter functions into a single function 

"""
NEED DOCUMENTATION
"""
function filter_data(trace::Experiment{T}; kwargs...) where {T<:Real}
    data = deepcopy(trace)
    filter_data!(data; kwargs...)
    return data
end

function filter_data!(trace::Experiment{T}; 
        freq_start=1.0, freq_stop = 55.0, bandwidth = 10.0,
        mode = :Lowpass, method = :Chebyshev2, 
        pole=8, ripple = 15.0, attenuation = 100.0
    ) where {T<:Real}

    #Determine the filter response
    if mode == :Lowpass
        responsetype = Lowpass(freq_stop; fs=1 / trace.dt)
    elseif mode == :Highpass
        responsetype = Highpass(freq_start; fs=1 / trace.dt)
    elseif mode == :Bandpass
        responsetype = Bandpass(freq_start, freq_stop, fs = 1/trace.dt)
    elseif mode == :Bandstop
        responsetype = Bandstop(freq_start, freq_stop, fs = 1/trace.dt)
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

#=
struct FilterCasette{T}
    filterMode::Vector{Symbol}
    filterMethods::Vector{Symbol}
    filterRNG::Vector{Tuple{T,T}}
end

function filter_data!(trace::Experiment{T}, casette)

end
=#
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

"""
Return CWT returns a 4D CWT array

cwt = [swp, data, wavelets, chs]

Will always return cwt
"""
function cwt_filter!(trace::Experiment{T}; wave=cDb2, β::Int64=2, 
    period_window::Tuple{Int64,Int64} = (1, 9), 
    level_window::Tuple{T, T} = (-Inf, Inf),
) where T <: Real
    c = wavelet(wave, β=β)
    n_wavelets = 100 #This will be changed posteriorly
    cwt_wave = zeros(size(trace,1), size(trace, 2), n_wavelets, size(trace,3))
    for swp = 1:size(trace, 1), ch = 1:size(trace, 3)
        y = ContinuousWavelets.cwt(trace[swp, :, ch], c)
        #It seems the only way to change 
        cwt_wave[swp,:, 1:size(y,2), ch] .= y
        reconstruct = zeros(size(y))
        if !any(period_window .== -1)
            reconstruct[:, period_window[1]:period_window[2]] .= y[:, period_window[1]:period_window[2]]
            trace.data_array[swp, :, ch] .= ContinuousWavelets.icwt(reconstruct, c, PenroseDelta()) |> vec
        else
            #zero all numbers outside of the level window
            outside_window = findall(level_window[1] .< y .< level_window[2])
            reconstruct[:, period_window[1]:period_window[2]] .= y[:, period_window[1]:period_window[2]]
            trace.data_array[swp, :, ch] .= ContinuousWavelets.icwt(reconstruct, c, PenroseDelta()) |> vec
        end
        n_wavelets = size(y, 2) #This allows us to use the posterior knowledge to change the array
    end
    return cwt_wave[:, :, (1:n_wavelets), :]
end

"""

"""
function dwt_filter!(trace::Experiment; wave=WT.db4, period_window::Tuple{Int64,Int64}=(1, 8))
    #In this case we have to limit the analyis to the window of dyadic time
    #This means that we can only analyze sizes if they are equal to 2^dyadic
    dyad_n = trunc(Int64, log(2, size(trace, 2)))
    dwt_wave = zeros(size(trace, 1), 2^dyad_n, size(trace, 3))
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
        println("Here")
        xt = dwt(x, wavelet(wave), dyad_n)
        dwt_wave[swp, :, ch] .= xt
        reconstruct = zeros(size(xt))
        reconstruct[2^period_window[1]:2^(period_window[2])] .= xt[2^period_window[1]:2^(period_window[2])]
        trace.data_array[swp, 1:2^dyad_n, ch] .= idwt(reconstruct, wavelet(wave))
    end
    return dwt_wave
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

function normalize(trace::Experiment; rng=(-1, 0))
    data = deepcopy(trace)
    normalize!(data)
    return data
end

function normalize_channel!(trace::ePhys.Experiment; rng=(-1, 0))
     if rng[1] < 0.0
          mins = minimum(minimum(data, dims = 2), dims = 1)
          trace.data_array ./= -mins
     else
          mins = maximum(maximum(data, dims = 2), dims = 1)
          trace.data_array ./= mins
     end
end

function normalize_channel(trace::ePhys.Experiment; rng=(-1, 0))
    data = deepcopy(trace)
    normalize_channel!(data)
    return data
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