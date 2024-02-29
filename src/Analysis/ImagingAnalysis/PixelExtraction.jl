collapse_dims(a) = dropdims(a, dims = (findall(size(a) .== 1)...,))

function zProject(exp::Experiment{TWO_PHOTON, T}) where T <: Real 
    px, py = exp.HeaderDict["framesize"]
    z_proj = (sum(exp, dims = 2) / size(exp, 2))
    return reshape(z_proj, (px, py))
end
 
frameAverage(data::Experiment{TWO_PHOTON, T}) where T<: Real = sum(exp, dims = 1)/size(exp,1)

function binarize(img, threshold = 0.5)
    return img .> threshold
end

function findROIcentroid(data, ROI_id::Int64)
    ROI_mask = getROImask(data, ROI_id)
    ROI_mask_coords = map(coord -> [coord[1],coord[2]], findall(ROI_mask.!= 0))
    ROI_coords_arr = hcat(ROI_mask_coords...)
    ROI_centroid = Tuple(mean(ROI_coords_arr, dims =2)')
    #Find the centroid by taking the mean
    return ROI_centroid
end

function findROIcentroid(data)
    ROI_mask = getROImask(data)
    ROI_centroids = map(i -> findROIcentroid(data, i), 1:maximum(ROI_mask))
    return ROI_centroids
end