collapse_dims(a) = dropdims(a, dims = (findall(size(a) .== 1)...,))

function zProject(exp::Experiment{TWO_PHOTON, T}) where T <: Real 
     px, py = exp.HeaderDict["framesize"]
     z_proj = (sum(exp, dims = 2) / size(exp, 2))
     return reshape(z_proj, (px, py))
 end
 
 frameAverage(data::Experiment{TWO_PHOTON, T}) where T<: Real = sum(exp, dims = 1)/size(exp,1)

 export zProject, frameAverage