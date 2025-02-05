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

"""
    baseline_als(y::Vector{Float64}; lam=1e7, p=0.25, niter=20) -> Vector{Float64}

Applies asymmetric least squares smoothing for baseline correction.
"""
function baseline_als(y::Vector{Float64}; lam=1e7, p=0.25, niter=20)
    L = length(y)
    n = L - 2
    rows, cols, vals = Int[], Int[], Float64[]

    for i in 1:n
        push!(rows, i); push!(cols, i); push!(vals, 1.0)
        push!(rows, i); push!(cols, i + 1); push!(vals, -2.0)
        push!(rows, i); push!(cols, i + 2); push!(vals, 1.0)
    end

    D = sparse(rows, cols, vals, n, L)
    DTD = lam * (D' * D)
    w = ones(Float64, L)
    z = copy(y)
    
    for _ in 1:niter
        W = spdiagm(0 => w)
        Z = W + DTD
        z = Z \ (w .* y)
        w = p .* (y .> z) .+ (1 - p) .* (y .< z)
    end

    return z
end

"""
    roi_function_fit(baselined_trace::Vector{Float64}, stim_frame::Int; starting_frame::Int=stim_frame) -> Tuple{Vector{Float64}, Float64, Float64, Float64, Float64}

Fits a function to the baselined trace.
"""
function roi_function_fit(baselined_trace::Vector{Float64}, stim_frame::Int; starting_frame::Int = stim_frame)
    window_start = starting_frame - 5
    x_window = ((window_start + 1):length(baselined_trace)) .- window_start
    x_window = x_window ./ 8.33
    y_window = baselined_trace[window_start + 1:end]

    lower_bounds = [-1.0, 0.001, 0.001, -1.0]
    upper_bounds = [5.0, 200.0, 100.0, 1.0]
    p0 = [(lb + ub) / 2 for (lb, ub) in zip(lower_bounds, upper_bounds)]
    model(x, p) = p[1] .* (1 .- exp.(-x ./ p[2])) .* exp.(-x ./ p[3]) .+ p[4]

    fit_result = try
        curve_fit(model, x_window, y_window, p0, lower = lower_bounds, upper = upper_bounds)
    catch
        nothing
    end

    optimal = isnothing(fit_result) ? [5.0, 100.0, 100.0, 1.0] : fit_result.param
    fit_trace = map(x -> max(x, 0.0), model(x_window, optimal))

    return fit_trace, optimal[1], optimal[2], optimal[3], optimal[4]
end

"""
    process_file(tiff_stack_path::String, roi_size::Int; stim_side::String="left", stim_frame::Int=200, ma_tail_start::Int=485, loading_bar::Bool=true)

Processes the TIFF stack by splitting it into ROIs and applying processing.
"""
function process_file(tiff_stack_path::String, roi_size::Int; stim_side::String = "left", stim_frame::Int = 200, ma_tail_start::Int = 485, loading_bar::Bool = true)
    # Load image using readImage (which handles normalization)
    raw_stack = readImage(tiff_stack_path)

    # Convert to Float64 if needed
    stack = Float64.(raw_stack)

    # Determine pixel splits
    image_size = size(stack)[2:3]
    xsplits, ysplits = pixel_splits(image_size, roi_size)

    roi_num_counter = 0
    n_x = length(xsplits) - 1
    n_y = length(ysplits) - 1
    total = n_x * n_y
    p = loading_bar ? Progress(total, 1) : nothing

    for i in 1:n_x
        for j in 1:n_y
            roi_array = stack[:, xsplits[i] + 1:xsplits[i + 1], ysplits[j] + 1:ysplits[j + 1]]
            pooled_trace = [mean(view(roi_array, k, :, :)) for k in 1:size(roi_array, 1)]
            baselined_trace = baseline_als(pooled_trace)
            fit_trace, a, tau1, tau2, c = roi_function_fit(baselined_trace, stim_frame)
            # Here you can store or process the fit_trace and parameters as needed
            roi_num_counter += 1
            if loading_bar
                next!(p)
            end
        end
    end
end