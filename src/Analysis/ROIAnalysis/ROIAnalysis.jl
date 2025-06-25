function findROIcentroids(array::Array{Int,2}, xlim::AbstractRange, ylim::AbstractRange)
    # Ensure the array is square
    size_x, size_y = size(array)
    @assert size_x == size_y "The input array must be square."

    # Calculate the number of elements along each axis
    n = size_x

    # Calculate the physical spacing between pixels
    dx = (xlim[2] - xlim[1])
    dy = (ylim[2] - ylim[1])

    # Find unique labels excluding 0 (assumed to be background)
    labels = setdiff(unique(array), [0])

    # Dictionary to store centroids: label => (x_centroid, y_centroid)
    centroids = Dict{Int, Tuple{Float64, Float64}}()

    for label in labels
        # Find indices where array equals the current label
        inds = findall(array .== label)

        if isempty(inds)
            continue  # Skip if no indices found for the label
        end

        # Initialize sums for x and y positions
        sum_x = 0.0
        sum_y = 0.0

        for I in inds
            i, j = Tuple(I)

            # Map array indices to physical coordinates
            x_pos = xlim[1] + (j - 0.5) * dx
            y_pos = ylim[1] + (i - 0.5) * dy

            sum_x += x_pos
            sum_y += y_pos
        end

        # Compute the centroid by averaging positions
        num_points = length(inds)
        centroid_x = sum_x / num_points
        centroid_y = sum_y / num_points

        # Store the centroid in the dictionary
        centroids[label] = (centroid_x, centroid_y)
    end

    return centroids
end

findROIcentroids(exp) = findROIcentroids(getROImask(exp), exp.HeaderDict["xrng"], exp.HeaderDict["yrng"])