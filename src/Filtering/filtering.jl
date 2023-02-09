#TODO: Eventually I want to add all filter functions into a single function 

"""
    Filters data in the `exp` object using a digital filter.

    # Parameters
    - exp : Experiment{T} 
        an object containing the data to be filtered
    - freq_start : float (optional) 
        start frequency for the filter, default is 1.0
    - freq_stop : float (optional) 
        stop frequency for the filter, default is 55.0
    - bandwidth : float (optional) 
        bandwidth for the filter, default is 10.0
    - mode : Symbol (optional) 
        filter mode, can be :Lowpass, :Highpass, :Bandpass, :Bandstop, default is :Lowpass
    - method : Symbol (optional) 
        method used to design the filter, can be :Butterworth, :Chebyshev1, :Chebyshev2, :Elliptic, default is :Chebyshev2
    - pole : int (optional) 
        number of poles for the filter, default is 8
    - ripple : float (optional) 
        ripple for the filter, default is 15.0
    - attenuation : float (optional) 
        attenuation for the filter, default is 100.0
    - filter_channels : int or Vector{String} (optional) 
        channels to be filtered, can be an integer index or a vector of channel names, default is -1 (all channels)
    # Returns
    - The input exp object is modified in place.

    # Example
    ```julia
    exp = Experiment(data_array, dt)
    filter_data!(exp, freq_start=5, freq_stop=10, mode=:Bandpass)
    ```
"""
function filter_data(exp::Experiment{T}; kwargs...) where {T<:Real}
    data = deepcopy(exp)
    filter_data!(data; kwargs...)
    return data
end

function filter_data!(exp::Experiment{T};
    freq_start=1.0, freq_stop=55.0, bandwidth=10.0,
    mode=:Lowpass, method=:Chebyshev2,
    pole=4, ripple=50.0, attenuation=100.0,
    filter_channels=-1
) where {T<:Real}

    #Determine the filter response
    if mode == :Lowpass
        responsetype = Lowpass(freq_stop; fs=1 / exp.dt)
    elseif mode == :Highpass
        responsetype = Highpass(freq_start; fs=1 / exp.dt)
    elseif mode == :Bandpass
        responsetype = Bandpass(freq_start, freq_stop, fs=1 / exp.dt)
    elseif mode == :Bandstop
        responsetype = Bandstop(freq_start, freq_stop, fs=1 / exp.dt)
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
        digital_filter = iirnotch(freq, bandwidth, fs=1 / exp.dt)
    else
        digital_filter = digitalfilter(responsetype, designmethod)
    end

    if filter_channels == -1
        filter_channels = 1:size(exp, 3)
    elseif isa(filter_channels, Vector{String}) #Do this for channel names input 
        filter_channels = findall(filter_channels .== data.chNames)
    end
    for swp in axes(exp, 1)
        for ch in axes(exp, 3)
            if ch in filter_channels
                exp.data_array[swp, :, ch] .= filt(digital_filter, exp[swp, :, ch])
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

function filter_data!(exp::Experiment{T}, casette)

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
    if t_post = 4.0, wave = Morlet(0.50π)
    if t_post = 2.0, wave = Morlet(0.25π)
    if t_post = 1.0, wave = Morlet(0.25π)
"""
function cwt_filter!(exp::Experiment{T};
    wave=Morlet(1.0π), β=2.0,
    period_window::Tuple{Real,Real}=(-1, -1),
    power_window::Tuple{T,T}=(0.0, 1.0),
    inverseStyle = NaiveDelta(),
    real_or_abs = :absolute,
    #Here are some other Wavelet options
    kwargs...
) where {T<:Real}
    c = wavelet(wave; β=β, kwargs...)
    n_wavelets = 200 #This will be changed posteriorly
    if isa(wave, Morlet)
        cwt_wave = zeros(ComplexF64, size(exp, 1), size(exp, 2), n_wavelets, size(exp, 3))
        #println("Complex")
    else
        cwt_wave = zeros(size(exp, 1), size(exp, 2), n_wavelets, size(exp, 3))
    end

    for swp = axes(exp, 1), ch = axes(exp, 3)
        y = ContinuousWavelets.cwt(exp[swp, :, ch], c)
        #It seems the only way to change
        if eltype(y) == ComplexF64
            if real_or_abs == :absolute
                check_y = abs.(y)
                println("Absolute")
            elseif real_or_abs == :real
                check_y = real.(y)
                println("Real")
            end
        else
            check_y = y
        end
        #println(minimum(check_y))
        #println(maximum(check_y))
        #we can fill the reconstruct with the lowest number 
        #println(min_val)
        min_val = y[argmin(abs.(y))]
        reconstruct = fill(min_val, size(y))
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
        reconstruct = reconstruct .* outside_window
        cwt_wave[swp, :, axes(reconstruct, 2), ch] .= reconstruct
        exp.data_array[swp, :, ch] .= ContinuousWavelets.icwt(reconstruct, c, inverseStyle) |> vec
        n_wavelets = size(reconstruct, 2) #This allows us to use the posterior knowledge to change the array
    end
    return cwt_wave[:, :, (1:n_wavelets), :]
end

function cwt_filter(exp::Experiment{T}; kwargs...) where {T<:Real}
    data = deepcopy(exp)
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
function dwt_filter!(exp::Experiment; wave=WT.db4, period_window::Tuple{Int64,Int64}=(1, 8))
    #In this case we have to limit the analyis to the window of dyadic time
    #This means that we can only analyze sizes if they are equal to 2^dyadic
    dyad_n = trunc(Int64, log(2, size(exp, 2)))
    dwt_wave = zeros(size(exp, 1), 2^dyad_n, size(exp, 3))
    #println(2^dyad_n)
    if (length(exp.t) - 2^dyad_n) > 0
        @warn "Sampling rate is non-dyadic. Will result in $(length(exp.t) - 2^dyad_n) leftover points"
    end
    if period_window[2] > dyad_n
        println("Period Window larger than dyad")
        period_window = (period_window[1], dyad_n)
    end
    for swp = axes(exp, 1), ch = axes(exp, 3)
        x = exp[swp, 1:2^dyad_n, ch]
        println("Here")
        xt = dwt(x, wavelet(wave), dyad_n)
        dwt_wave[swp, :, ch] .= xt
        reconstruct = zeros(size(xt))
        reconstruct[2^period_window[1]:2^(period_window[2])] .= xt[2^period_window[1]:2^(period_window[2])]
        exp.data_array[swp, 1:2^dyad_n, ch] .= idwt(reconstruct, wavelet(wave))
    end
    return dwt_wave
end

function dwt_filter(exp::Experiment; kwargs...)
    #In this case we have to limit the analyis to the window of dyadic time
    #This means that we can only analyze sizes if they are equal to 2^dyadic
    #We can fix this by taking a 
    data = deepcopy(exp)
    dwt_filter!(data; kwargs...)
    return data
end

"""
This is from the adaptive line interface filter in the Clampfit manual

This takes notch filters at every harmonic

#Stimulus artifacts have a very specific harmonic
250, 500, 750, 1000 ... 250n
"""
function EI_filter(exp; reference_filter=60.0, bandpass=10.0, cycles=5)
    data = deepcopy(exp)
    for cycle in 1:cycles
        band!(data, center=reference_filter * cycle, std=bandpass)
    end
    return data
end

function EI_filter!(exp; reference_filter=60.0, bandpass=10.0, cycles=5)
    for cycle in 1:cycles
        notch_filter!(exp, center=reference_filter * cycle, std=bandpass)
    end
end

function normalize!(exp::Experiment; rng=(-1, 0), normalize_by=:channel)
    for swp in axes(exp, 1)
        for ch in axes(exp, 3)
            if rng[1] < 0
                exp.data_array[swp, :, ch] .= (exp[swp, :, ch] ./ minimum(exp[swp, :, ch], dims=2))
            else
                exp.data_array[swp, :, ch] .= (exp[swp, :, ch] ./ maximum(exp[swp, :, ch], dims=2))
            end
        end
    end
end

function normalize(exp::Experiment; rng=(-1, 0), dims=2)
    data = deepcopy(exp)
    normalize!(data)
    return data
end

function normalize_channel!(exp::ePhys.Experiment; rng=(-1, 0))
    if rng[1] < 0.0
        mins = minimum(minimum(data, dims=2), dims=1)
        exp.data_array ./= -mins
    else
        mins = maximum(maximum(data, dims=2), dims=1)
        exp.data_array ./= mins
    end
end

function normalize_channel(exp::ePhys.Experiment; rng=(-1, 0))
    data = deepcopy(exp)
    normalize_channel!(data)
    return data
end

function rolling_mean(exp::Experiment; window::Int64=10)
    data = deepcopy(exp)
    for swp in axes(exp, 1), ch in axes(exp, 3)
        for i in 1:window:size(data, 2)-window
            data.data_array[swp, i, ch] = sum(data.data_array[swp, i:i+window, ch]) / window
        end
    end
    return data
end