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
        pole=8, ripple = 15.0, attenuation = 100.0, 
        filter_channels = -1
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

    if filter_channels == -1
        filter_channels = 1:size(trace,3)
    elseif isa(filter_channels, Vector{String}) #Do this for channel names input 
        filter_channels = findall(filter_channels  .== data.chNames)
    end
    for swp in 1:size(trace, 1)
        for ch in 1:size(trace, 3)
            if ch in filter_channels
                trace.data_array[swp, :, ch] .= filt(digital_filter, trace[swp, :, ch])
            end
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


"""
Return CWT returns a 4D CWT array

cwt = [swp, data, wavelets, chs]

Will always return cwt

Here are some useful settings for filtering artifacts
    β = 1.0
    wave = cDb2
    level_window = (0.5, 1.0)
    period_window = (20, 60)
"""
function cwt_filter!(trace::Experiment{T}; 
    wave=cDb2, β = 2.0, 
    period_window::Tuple{Int64,Int64} = (-1, -1), 
    power_window::Tuple{T, T} = (0.0, 1.0),
    inverseStyle = NaiveDelta(), 
    #Here are some other Wavelet options
    kwargs...
) where T <: Real
    c = wavelet(wave; β = β, kwargs...)
    n_wavelets = 200 #This will be changed posteriorly
    if isa(wave, Morlet)
        cwt_wave = zeros(ComplexF64, size(trace,1), size(trace, 2), n_wavelets, size(trace,3))
        #println("Complex")
    else
        cwt_wave = zeros(size(trace,1), size(trace, 2), n_wavelets, size(trace,3))
    end

    for swp = 1:size(trace, 1), ch = 1:size(trace, 3)
        y = ContinuousWavelets.cwt(trace[swp, :, ch], c)
        #It seems the only way to change
        if eltype(y)==ComplexF64
            check_y = real.(y) 
        else
            check_y = y
        end
        #println(minimum(check_y))
        #println(maximum(check_y))

        reconstruct = zeros(eltype(y), size(y))
        #println(eltype(y))
        if period_window[1] == -1 && period_window[2] != -1
            reconstruct[:, 1:period_window[2]] .= y[:, 1:period_window[2]]
        elseif period_window[1] != -1 && period_window[2] == -1
            reconstruct[:, period_window[1]:end] .= y[:, period_window[1]:end]
        elseif period_window[1] != -1 && period_window[2] != -1 
            reconstruct[:, period_window[1]:period_window[2]] .= y[:, period_window[1]:period_window[2]]
        else
            reconstruct = y
        end
        
        #standardize the check
        norm_y = standardize(UnitRangeTransform, check_y)
        outside_window = (power_window[1] .<= norm_y .<= power_window[2])
        reconstruct = reconstruct.*outside_window
        cwt_wave[swp,:, 1:size(reconstruct,2), ch] .= reconstruct
        trace.data_array[swp, :, ch] .= ContinuousWavelets.icwt(reconstruct, c, inverseStyle) |> vec
        n_wavelets = size(reconstruct, 2) #This allows us to use the posterior knowledge to change the array
    end
    return cwt_wave[:, :, (1:n_wavelets), :]
end

function cwt_filter(trace::Experiment{T}; kwargs...) where T <: Real
    data = deepcopy(trace)
    cwt = cwt_filter!(data; kwargs...)
    return data, cwt
end

# Wavelet utility funtcions 
function CWTprocess(y)
    y2 = log.(2, (abs.(y) .^ 2))
    #return standardize(UnitRangeTransform, y2)'
    return y2'
end

function coi_factor(i, j, n)
    # Compute the normalized distance from the center of the signal
    r = (i - 1) / (n / 2)
    # Compute the cone of influence factor
    return 1.0 - abs(r)
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

function normalize!(trace::Experiment; rng=(-1, 0), normalize_by = :channel)
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

function normalize(trace::Experiment; rng=(-1, 0), dims = 2)
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