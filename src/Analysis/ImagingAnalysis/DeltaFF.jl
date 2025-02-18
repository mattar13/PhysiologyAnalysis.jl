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
function baseline_als(y::Vector{T}; lam::T=1e7, p::T=0.075, niter::Int=20) where T<:Real
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
        w = [ y[i] > z[i] ? p : (y[i] < z[i] ? 1-p : 0.5 ) for i in 1:L ]
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
Perform overall baseline correction using ALS and centered Moving Average.

- `stim_frame=0`: Disables stimulus-based normalization.
- `ma_tail_start=0`: Disables ignored region in the moving average correction.
- `window`: Size of the moving average window.

Returns a baseline-corrected dF/F trace.
"""
function baseline_trace(trace::AbstractVector{T}; 
    window::Int=15, #This is the window of the moving average for dF
    kwargs...
) where T<:Real

    # Normalize using pre-stimulus baseline if stim_frame is set
    baseline_divisor = mean(trace)  # Default to global mean if no stimulus

    F0 = trace ./ baseline_divisor

    # Apply Asymmetric Least Squares (ALS) smoothing
    drift = baseline_als(F0; kwargs...)
    dF = F0 .- drift

    dFoF = moving_average(dF; window = window)

    # Compute final dF/F trace
    #dFoF = baselined .- ma_base
    return dFoF
end


# 2) Stack functions ============================================================#
"""
    moving_average(A, window)

Compute a centered moving average along each column of a 2D array `A`
using ImageFiltering.jl's `imfilter`. The kernel is a uniform averaging filter
of length `window`. A symmetric border condition (`Reflect()`) is used so that the
filter is applied in a centered manner. (For best centering, use an odd window size.)

Returns an array of the same size as `A`.
"""

function moving_average(A::AbstractMatrix{T}; window::Int = 15, border = "reflect") where T<:Real
    # Create a 1D averaging kernel (as a column vector)
    kernel = reshape(ones(T, window) ./ window, (window, 1))
    B = imfilter(T, A, kernel,border)
    return B  # or return B_offset
end

"""
    baseline_stack(stack; kwargs...)
"""
function baseline_stack(stack::Matrix{T}; window = 15, kwargs...) where T<:Real
    #First I want to calculate the global drift
    global_trace = mean(stack, dims = 1)[1,:]
    baseline_divisor = mean(global_trace) # Default to global mean if no stimulus
    global_F0 = global_trace ./ baseline_divisor
    #Apply Asymmetric Least Squares (ALS) smoothing
    drift = baseline_als(global_F0; kwargs...)
    
    #Apply the moving average to the stack
    F0_stack = stack ./ baseline_divisor #Adjust the stack by the global mean
    dF = F0_stack .- reshape(drift, (1, length(drift))) #calculate the delta F0
    dFoF = moving_average(dF; window = window)
    return dFoF
end