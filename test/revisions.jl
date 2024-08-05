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

#%% ╔═╡ I want to use the peakfinder algorithim
pks, vals = findmaxima(delta_f_f_grn_trace, 40)
pks = pks[vals.<10.0]
vals = vals[vals.<10.0]

#Do the crosscorrelation analysis
using StatsBase

#%% ╔═╡ Plot the results of the peak finding
fig1 = Figure(size = (1000, 800))
ax1 = GLMakie.Axis(fig1[1,1:2], title = "Green Trace")#, aspect = 1.0)
ax2 = GLMakie.Axis(fig1[2,1:2], title = "Green Trace")#, aspect = 1.0)
lines!(ax1, data2P.t, delta_f_f_grn_trace, color = :green)
lines!(ax2, data2P.t, delta_f_f_red_trace, color = :red)

scatter!(ax1, data2P.t[pks], vals, color = :red, markersize = 10.0)

ax3a = GLMakie.Axis(fig1[3,1], title = "Section")#, aspect = 1.0)
ax3b = GLMakie.Axis(fig1[3,2], title = "Section")#, aspect = 1.0)

#Plot all of the wave events on top of each other (using the 40)
event_length = 50
new_t = (-event_length/2:event_length/2-1) * data2P.dt
red_sect_arr = zeros(event_length, length(pks))
for (i, pk) in enumerate(pks)
    t = pk#data2P.t[pk]
    if t-event_length/2 < 0
        t_start = 1
    else
        t_start = round(Int64, t-event_length/2)
    end

    if t+event_length/2 > length(delta_f_f_grn_trace)
       t_end = length(delta_f_f_grn_trace) 
    else
        t_end = round(Int64, t+event_length/2-1)
    end

    t_rng = t_start:t_end
    grn_sect = delta_f_f_grn_trace[t_rng]
    red_sect = delta_f_f_red_trace[t_rng]
    red_sect_arr[1:length(red_sect), i] = red_sect
    
    lines!(ax3a, new_t, grn_sect, color = :green)
    lines!(ax3b, new_t, red_sect, color = :red, alpha = 0.2)
end
red_sect_arr
red_mean = mean(red_sect_arr, dims = 2)[:,1]
lines!(ax3b, new_t, red_mean, color = :red)
fig1

#%% ╔═╡Plot the figure of the cells and their background
fig2 = Figure(size = (1000, 800))
ax2b = GLMakie.Axis(fig2[1,1], title = "Green Channel", aspect = 1.0)
ax2c = GLMakie.Axis(fig2[3,1], title = "Red Channel", aspect = 1.0)
ax3b = GLMakie.Axis(fig2[1,2:3], title = "Green Trace")#, aspect = 1.0)
ax3c = GLMakie.Axis(fig2[3,2:3], title = "Red Trace")#, aspect = 1.0)

ax3b.ylabel = "Delta f/f"
ax3c.ylabel = "Delta f/f"
ax3c.xlabel = "Time (ms)"
hm2a = heatmap!(ax2b, xlims, ylims, delta_f_f_grn_zstack[:,:,1], colormap = :viridis, colorrange = (0.0, maximum(delta_f_f_grn_trace)))
hm2b = heatmap!(ax2c, xlims, ylims, delta_f_f_red_zstack[:,:,1], colormap = :CMRmap, colorrange = (0.0, maximum(delta_f_f_red_trace)))
Colorbar(fig2[2, 1], hm2a, vertical = false, width = @lift Fixed($(pixelarea(ax2b.scene)).widths[1]))
Colorbar(fig2[4, 1], hm2b, vertical = false, width = @lift Fixed($(pixelarea(ax2c.scene)).widths[1]))
lines!(ax3b, data2P.t, delta_f_f_grn_trace, color = :green)
lines!(ax3c, data2P.t, delta_f_f_red_trace, color = :red)

fig1



#%% ╔═╡ Lets use CellPose to label cells
#Pkg.build("PyCall")
using PyCall, Conda
model = cellpose_model()
mask, flow, style, diam = model.eval(grn_zproj)
data2P.HeaderDict["ROIs"] .= vec(mask)

#%% ╔═╡Plot the figure
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


roi_mask = getROImask(data2P)
idx = 1
coords = findall(roi_mask .== idx)

for (k,v) in data2P.HeaderDict
    println(k)
end