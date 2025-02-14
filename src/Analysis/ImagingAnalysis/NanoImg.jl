#These functions are used to do some pixel extraction and ROI analysis
using DataFrames
using SparseArrays
"""
    pixel_splits(image_size::Tuple{Int, Int}, roi_size::Int) -> Tuple{Vector{Int}, Vector{Int}}

Determines the pixel splitting indices for the image based on `roi_size`.
"""
function pixel_splits(image_size::Tuple{Int, Int}, roi_size::Int)
    x_pixels, y_pixels = image_size
    xs = collect(0:roi_size:x_pixels)
    ys = collect(0:roi_size:y_pixels)

    if xs[end] < x_pixels
        push!(xs, x_pixels)
    end
    if ys[end] < y_pixels
        push!(ys, y_pixels)
    end

    return (xs, ys)
end

#############################
# 1. Baseline Correction
#############################

"""
    baseline_als(y; lam=1e7, p=0.25, niter=20)

Asymmetric Least Squares baseline correction.

Uses a second–difference regularization (similar to Eilers & Boelens, 2005).
"""
function baseline_als(y::Vector{T}; lam::T=1e7, p::T=0.25, niter::Int=20) where T<:Real
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
function moving_average(a::Vector{T}, window::Int) where T<:Real
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
Compute a centered moving–average baseline.

If `early_frame=0` and `late_frame=0`, applies a moving average across the entire trace.
Otherwise, interpolates linearly between `early_frame` and `late_frame`.

- `window`: Size of the moving average window.
- `early_frame`: Start of the ignored stimulation window (0 to disable).
- `late_frame`: End of the ignored stimulation window (0 to disable).
"""
function baseline_ma(trace::Vector{T}; 
    window::Int=40, 
    early_frame::Int=0, 
    late_frame::Int=0) where T<:Real
    L = length(trace)
    ma_bl = moving_average(trace, window)  # Compute centered moving average

    # If early_frame and late_frame are set (not 0), interpolate linearly
    if early_frame > 0 && late_frame > 0
        linear_fill = range(ma_bl[early_frame], stop=ma_bl[late_frame], length=late_frame-early_frame+1)
        ma_bl[early_frame:late_frame] = collect(linear_fill)
    end

    return ma_bl
end

"""
Perform overall baseline correction using ALS and centered Moving Average.

- `stim_frame=0`: Disables stimulus-based normalization.
- `ma_tail_start=0`: Disables ignored region in the moving average correction.
- `window`: Size of the moving average window.

Returns a baseline-corrected dF/F trace.
"""
function baseline_trace(trace::Vector{T}; 
    stim_frame::Int=0, 
    ma_tail_start::Int=0, 
    window::Int=40) where T<:Real

    L = length(trace)

    # Normalize using pre-stimulus baseline if stim_frame is set
    if stim_frame > 30  # Ensure valid range
        baseline_divisor = mean(trace[stim_frame-30:stim_frame-5])
    else
        baseline_divisor = mean(trace)  # Default to global mean if no stimulus
    end

    trace_adj = trace ./ baseline_divisor

    # Apply Asymmetric Least Squares (ALS) smoothing
    als_base = baseline_als(trace_adj)
    baselined = trace_adj .- als_base

    # Ignore ignored frames if stim_frame or ma_tail_start are 0
    early_frame = stim_frame > 0 ? stim_frame - 10 : 0
    late_frame = ma_tail_start > 0 ? ma_tail_start : 0

    ma_base = baseline_ma(baselined; window=window, 
        early_frame=early_frame, 
        late_frame=late_frame)

    # Compute final dF/F trace
    dFoF = baselined .- ma_base
    return dFoF
end