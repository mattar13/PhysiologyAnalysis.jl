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

    return dFoF
end