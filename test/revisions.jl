using Revise
using ElectroPhysiology, PhysiologyAnalysis
using Pkg; Pkg.activate("test")

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
using GLMakie, PhysiologyPlotting
data2P_fn = raw"D:\Data\Calcium Imaging\2024_07_24_OPN4g_P9\ca_img4004.tif"

#╔═╡Extract the image
data2P = readImage(data2P_fn);
xlims = data2P.HeaderDict["xrng"]
ylims = data2P.HeaderDict["yrng"]

deinterleave!(data2P) #This seperates the movies into two seperate movies

# ╔═╡Seperate the red and green channels
img_arr = get_all_frames(data2P)
red_zstack = img_arr[:,:,:,2]*100
grn_zstack = img_arr[:,:,:,1]*100
red_zproj = project(data2P, dims = (3))[:,:,1,2]
grn_zproj = project(data2P, dims = (3))[:,:,1,1]
red_trace = project(data2P, dims = (1,2))[1,1,:,2]
grn_trace = project(data2P, dims = (1,2))[1,1,:,1]
delta_f_f_red_zstack = deltaF_F(red_zstack; voxel_z = 200, mode = "symmetric")
delta_f_f_grn_zstack = deltaF_F(grn_zstack; voxel_z = 200, mode = "symmetric")
delta_f_f_red_trace = mean(delta_f_f_red_zstack, dims = (1,2))[1,1,:]
delta_f_f_grn_trace = mean(delta_f_f_grn_zstack, dims = (1,2))[1,1,:]

#%% ╔═╡ Lets use CellPose to label cells
#Pkg.build("PyCall")
using PyCall, Conda
model = cellpose_model()
mask, flow, style, diam = model.eval(grn_zproj)
data2P.HeaderDict["ROIs"] .= vec(mask)

# Can we find out the centroid of each img from ROIs?
mask = getROImask(data2P)
centroids = findROIcentroids(data2P)

# ╔═╡Plot the figure
fig2 = Figure(size = (1200, 500))
ax1 = GLMakie.Axis(fig2[1,1], title = "Green Channel", aspect = 1.0)
ax2 = GLMakie.Axis(fig2[1,2], title = "Fluor")

hm1 = heatmap!(ax1, xlims, ylims, grn_zproj, colormap = :viridis, colorrange = (0.0, maximum(grn_zproj)/25))
hm2 = contour!(ax1, xlims, ylims, mask, color = :red)
lines!(ax2, data2P.t, grn_trace)

for idx in 1:maximum(data2P.HeaderDict["ROIs"])
    println(idx)
    ROI_data = data2P[getROIindexes(data2P, idx), :, :]
    ROI_trace = mean(ROI_data, dims = 1)[1,:,1]
    lines!(ax2, data2P.t, ROI_trace)
end

for (k, v) in centroids
    println("$k - $v")
    scatter!(ax1, (v[2], v[1]))
end
fig2

#%%
for (k,v) in data2P.HeaderDict
    println(k)
end