 #############################
# 1. Baseline Correction
#############################

"""
    baseline_als(y; lam=1e7, p=0.25, niter=20)

Asymmetric Least Squares baseline correction.
# Parameters
- `y::Vector{T}`: The input signal to be corrected, where `T` is a subtype of `Real`.
- `lam::T=1e7`: The smoothing parameter (λ) that controls the smoothness of the baseline. A higher value enforces a smoother baseline by penalizing large differences between consecutive points. Typical values range from 10⁴ to 10⁹. Adjust `lam` based on the desired smoothness; larger values yield smoother baselines.
- `p::T=0.25`: The asymmetry parameter that determines the weighting between positive and negative residuals. It should be in the range (0, 1). For signals with positive peaks, `p` is usually set between 0.001 and 0.1. A smaller `p` reduces the influence of positive deviations, making the baseline less sensitive to peaks.
- `niter::Int=20`: The number of iterations for the algorithm to perform. Each iteration refines the baseline estimate. Typically, `niter` is set between 10 and 20. More iterations can lead to a more accurate baseline but will increase computation time.

If the peaks are very sudden, p gives you the best chance of fixing it
Uses a second–difference regularization (similar to Eilers & Boelens, 2005).
"""
function baseline_als(y::Vector{T}; 
    lam::T=1e7, assym::T=0.25, niter::Int=20
) where T<:Real
    L = length(y)
    # Construct D as an L×(L–2) difference operator so that
    # (D*y)[i] = y[i] - 2y[i+1] + y[i+2]
    D = spzeros(T, L, L-2)
    for j in 1:(L-2)
        D[j, j] = 1
        D[j+1, j] = -2
        D[j+2, j] = 1
    end
    # Compute penalty matrix (L×L)
    Dmat = lam * (D * transpose(D))
    
    w = ones(T, L)
    z = similar(y)
    for _ in 1:niter
        # Create diagonal weight matrix
        W = spdiagm(0 => w)
        Z = W + Dmat
        z .= Z \ (w .* y)
        # Update weights elementwise
        w = [ y[i] > z[i] ? assym : (y[i] < z[i] ? 1-assym : 0.5 ) for i in 1:L ]
    end
    return z
end


"""
Compute a centered moving average with a given window size.

- `window`: Number of points to average over.
- Returns a vector of the same length as `a`, with the moving average centered.
"""
function moving_average(a::Vector{T}; window::Int = 15) where T<:Real
    L = length(a)
    half_window = div(window, 2)  # Centering the window

    # Compute forward moving sum (standard cumulative sum trick)
    cs = cumsum(a)
    mavg = zeros(T, L)

    # Compute moving average in a centered manner
    for i in 1:L
        left = max(i - half_window, 1)
        right = min(i + half_window, L)
        mavg[i] = sum(a[left:right]) / (right - left + 1)
    end

    return mavg
end

"""
Compute a centered moving median with a given window size.

- `window`: Number of points to use for the median calculation.
- Returns a vector of the same length as `a`, with the moving median centered.
"""
function median_filter(a::Vector{T}; window::Int = 15) where T<:Real
    L = length(a)
    half_window = div(window, 2)  # Centering the window
    result = zeros(T, L)

    # Compute moving median in a centered manner
    for i in 1:L
        left = max(i - half_window, 1)
        right = min(i + half_window, L)
        result[i] = median(a[left:right])
    end

    return result
end

"""
    linear_fill!(a, start, stop)

Linearly interpolate the entries of `a` from index `start` to `stop`, inclusive,
so that `a[start]` and `a[stop]` remain fixed and all in-between points lie 
on the straight line between them.

This version:
  • works on any AbstractVector of real numbers,
  • uses a `LinRange` under the hood for clarity,
  • only writes the slice once.
"""
function linear_fill(a::AbstractVector{<:Real}, start::Integer, stop::Integer)
    @assert 1 ≤ start < stop ≤ length(a) "start and stop must satisfy 1 ≤ start < stop ≤ length(a)"
    # compute a range of exactly (stop-start+1) points
    vals = LinRange(a[start], a[stop], stop - start + 1)
    return vals
end

function linear_fill!(a::AbstractVector{<:Real}, start::Integer, stop::Integer)
    @assert 1 ≤ start < stop ≤ length(a) "start and stop must satisfy 1 ≤ start < stop ≤ length(a)"
    # compute a range of exactly (stop-start+1) points
    vals = LinRange(a[start], a[stop], stop - start + 1)
    a[start:stop] .= vals
    return a
end

"""
Perform overall baseline correction using ALS and centered Moving Average.

- `stim_frame=0`: Disables stimulus-based normalization.
- `ma_tail_start=0`: Disables ignored region in the moving average correction.
- `window`: Size of the moving average window.

Returns a baseline-corrected dF/F trace.
"""
function baseline_trace(trace::AbstractVector{T}; 
    window::Int = 20,
    stim_frame = nothing,
    baseline_divisor_start = 20, baseline_divisor_end = 5,
    linear_fill_start = nothing, linear_fill_end = nothing,
    warning = false,
    kwargs...
) where T<:Real

    # Normalize using pre-stimulus baseline if stim_frame is set
    if !isnothing(stim_frame) #We want to normalize using the pre-stimulus baseline
        if isnothing(baseline_divisor_end)
            baseline_divisor_end = length(trace)    
        end
        # Calculate the mean of the trace before the stimulus frame
        @assert baseline_divisor_start > baseline_divisor_end "baseline_divisor_start must be greater than baseline_divisor_end"
        baseline_divisor_start_idx = max(1, stim_frame - baseline_divisor_start)
        if stim_frame - baseline_divisor_start < 1 && warning
            @warn "Baseline Divisor Start exceeds the Stimulus frame"
            println("\t stim_frame: $(stim_frame), \n\t baseline_divisor_start: $(baseline_divisor_start)")
        end

        baseline_divisor_end_idx = max(1, stim_frame - baseline_divisor_end)
        if stim_frame - baseline_divisor_end < 1 && warning
            @warn "Baseline Divisor End exceeds the Stimulus frame"
            println("\t stim_frame: $(stim_frame), \n\t baseline_divisor_end: $(baseline_divisor_end)")
        end
        baseline_divisor = mean(trace[baseline_divisor_start_idx:baseline_divisor_end_idx])
    else #We want to normalize using the entire trace
        baseline_divisor = mean(trace)
    end
    F0 = trace ./ baseline_divisor
    
    # Apply Asymmetric Least Squares (ALS) smoothing
    drift = baseline_als(F0; kwargs...)
    baselined_trace = F0 .- drift

    #dFoF = moving_average(baselined_trace; window = window)
    #Enter in linear_fill in the sections where we need to
    ma = moving_average(baselined_trace; window = window)
    if !isnothing(linear_fill_start) && !isnothing(linear_fill_end) && !isnothing(stim_frame)
        linear_fill_start_idx = max(1, stim_frame - linear_fill_start)
        linear_fill_end_idx = min(length(trace), stim_frame + linear_fill_end)
        if linear_fill_start_idx == 1 && warning
            @warn "Linear fill starts at the beginning of the trace"
            println("\t linear_fill_start_idx: $linear_fill_start, \n\t stim_frame: $stim_frame")
        end
        if linear_fill_end_idx == length(trace) && warning
            @warn "Linear fill ends at the end of the trace"
            println("\t linear_fill_end_idx: $linear_fill_end, \n\t stim_frame: $stim_frame \n\t trace length: $(length(trace))")
        end
        linear_fill!(ma, linear_fill_start_idx, linear_fill_end_idx)
    end
    dFoF = baselined_trace - ma
    
    return dFoF
end

