using Revise
using ElectroPhysiology, PhysiologyAnalysis
using Pkg; Pkg.activate("test")

#%% ╔═╡This task is for filepath name and extraction using Dataframes.jl
example_fn = raw"D:\Data\Patching\2024_07_22_OPN4_P7\Cell5\24722019.abf"



#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
using GLMakie, PhysiologyPlotting
data2P_fn = raw"G:\Data\Calcium Imaging\2024_07_24_OPN4_P9\ca_img4004.tif"

#╔═╡Extract the image
data2P = readImage(data2P_fn);
xlims = data2P.HeaderDict["xrng"]
ylims = data2P.HeaderDict["yrng"]
deinterleave!(data2P) #This seperates the movies into two seperate movies

#╔═╡Seperate the red and green channels
img_arr = get_all_frames(data2P)

grn_trace = project(data2P, dims = (1,2))[1,1,:,1]
red_trace = project(data2P, dims = (1,2))[1,1,:,2]

grn_zproj = project(data2P, dims = (3))[:,:,1,1]
red_zproj = project(data2P, dims = (3))[:,:,1,2]

#Pkg.build("PyCall")
using PyCall, Conda
model = cellpose_model()
mask, flow, style, diam = model.eval(grn_zproj)
data2P.HeaderDict["ROIs"] .= vec(mask)

#%% ╔═╡Plot the figure
fig = Figure(size = (1200, 500))
ax1 = GLMakie.Axis(fig[1,1], title = "Green Channel", aspect = 1.0)
ax2 = GLMakie.Axis(fig[1,2], title = "Fluor")

hm1 = heatmap!(ax1, xlims, ylims, grn_zproj, colormap = :viridis, colorrange = (0.0, maximum(grn_zproj)/25))
hm2 = contour!(ax1, xlims, ylims, mask, color = :red)
grn_trace
time
lines!(ax2, data2P.t, grn_trace)

for idx in 1:maximum(data2P.HeaderDict["ROIs"])
    println(idx)
    ROI_data = data2P[getROIindexes(data2P, idx), :, :]
    ROI_trace = mean(ROI_data, dims = 1)[1,:,1]
    lines!(ax2, data2P.t, ROI_trace)
end


roi_mask = getROImask(data2P)
idx = 1
coords = findall(roi_mask .== idx)

for (k,v) in data2P.HeaderDict
    println(k)
end