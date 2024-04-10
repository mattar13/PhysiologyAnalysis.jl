using Revise
using ElectroPhysiology, PhysiologyAnalysis
ElectroPhysiology.__init__()
using Pkg; Pkg.activate("test")
using GLMakie
using FileIO, ImageView, Images
using Statistics
import GLMakie.Axis

#Root file for all data
data_root = raw"F:\Data"
ca_imaging_root = joinpath(data_root, "Calcium Images")
patching_root = joinpath(data_root, "Patching")
patch_fn = joinpath(patching_root, "2024_02_12_WT\\Cell1\\24212004.abf")
ca_img_fn = joinpath(ca_imaging_root, "2024_02_12_WT\\Cell1\\cell1001.tif")

# Open the data
dataIC = readABF(patch_fn); downsample!(dataIC, 100.0); #Go through these functions much more extensively
data2P = readImage(ca_img_fn; sampling_rate = 1.45);
n_frames = size(data2P, 2)
dt = data2P.dt
t = (collect(1:n_frames).-1) .* dt

mov = get_all_frames(data2P) #get all frames as a 3D array
fluo = mean(mov, dims = (1,2))[1,1,:]
zproj = mean(mov, dims = 3)[:,:,1] #Can we save this directly? 
px_x, px_y = data2P.HeaderDict["framesize"]
xlims = 1:px_x; 
ylims = 1:px_y
zlims = 1:size(data2P,2)
# Lets do some data analysis

#get the threshold
dataIC.data_array
mean(dataIC.data_array, dims = 2)

thresholds = calculate_threshold(dataIC, dims = 2)
tseries = dataIC.t
spike_array = dataIC.data_array .> thresholds

#%% Extract and plot the amacrine cell recordings
fig1 = Figure(size = (500, 700))
ax1 = Axis(fig1[1,1])
ax2 = Axis(fig1[2,1], ylabel = "Fluorescence")
ax3 = Axis(fig1[3,1], xlabel = "time (ms)", ylabel = "Voltage (mV)")
rowsize!(fig1.layout, 1, Relative(2/3))

data2P.t
heatmap!(ax1, xlims, ylims, zproj)
lines!(ax2, data2P.t, fluo)
lines!(ax3, dataIC.t, dataIC.data_array[1, :, 1])
display(fig1)
#%% Conduct the analysis





#%% Plot the current clamp with the 
patch_fn2 = joinpath(patching_root, "2024_02_12_WT\\Cell1\\24212006.abf")
data2 = readABF(patch_fn2)
downsample!(data2, 1000.0) #Go through these functions much more extensively
size(data2)
# Extract and plot the amacrine cell recordings
fig2 = Figure(size = (1000, 1000))
ax1 = Axis(fig2[1,1])
ax2 = Axis(fig2[2,1])
for tr in 1:size(data2, 1)
    lines!(ax1, data2.t, data2.data_array[tr, :, 1])
    lines!(ax2, data2.t, data2.data_array[tr, :, 2])
end
display(fig2)
data.chUnits
fluo = mean(mov, dims = (1,2))[1,1,:]
GLMakie.lines!(ax31, time, fluo)

for i in 1:maximum(getROImask(data))
    ROIi_mask = getROImask(data, i, reshape_arr = false)[:,1]
    coords = findall(ROIi_mask .!= 0.0)
    ROI_trace = mean(data.data_array[coords, :, 1], dims = 1)[1,:]
    GLMakie.lines!(ax32, time, ROI_trace)
end

ax41 = GLMakie.Axis(fig[4,1])
ax42 = GLMakie.Axis(fig[4,2])
scatter!(ax41, centroids, markersize = 20.0)
scatter!(ax42, centroids, markersize = 20.0)

record(fig, "ANIMATION.mp4", 2:size(data,2)) do i 
     println(i)
     hm1[3] = mov[:,:,i]
     hm2[3] = ROIall_mov[:,:,i]
end
display(fig)