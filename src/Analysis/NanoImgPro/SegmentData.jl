#These functions are used to do some pixel extraction and ROI analysis
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