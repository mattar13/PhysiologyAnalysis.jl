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