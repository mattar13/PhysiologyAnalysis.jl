#%% Load Packages _________________________________________________________________________#
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyAnalysis
Pkg.activate("test")
using GLMakie
import GLMakie.Axis
using FileIO, ImageView, Images
using Statistics

save_fn = "D:/Data/Analysis/2024_02_12_WTCell1_IC2P.mp4"

#%%=[Point to this for data from cell 1]=========================================================#
data2P_fn = "D:/Data/Calcium Imaging/2024_02_12_WT/Cell1/cell1001.tif"
dataIC_fn = "D:/Data/Patching/2024_02_12_WT/Cell1/24212004.abf"

#%%=[Point to this for data from cell 3]=========================================================#
data2P_fn = "D:/Data/Calcium Imaging/2024_03_11_VglutGC6/Cell3/freeRec001.tif"

#%%=[Open and prepare 2P data]=========================================================#
fps = 1.45
data2P = readImage(data2P_fn; sampling_rate = fps, chName = "525 nm");
px_x, px_y = data2P.HeaderDict["framesize"]
xlims = data2P.HeaderDict["xrng"]; #Extract the x domain
ylims = data2P.HeaderDict["yrng"]; #Extract the y domain

t = data2P.t #Extract the t domain
#If every other frame is a different channel uncomment line below
t = t[1:2:end]

mov = get_all_frames(data2P)[:,:,1:2:end] #get all frames as a 3D array
fluo = mean(mov, dims = (1,2))[1,1,:]
zproj = maximum(mov, dims = 3)[:,:,1] #Can we save this directly? 

# Open and prepare the data______________________________________________________________#
dataIC = readABF(dataIC_fn); 
downsample!(dataIC, 1000.0); #Go through these functions much more extensively

# Create the plot_______________________________________________________________________#
fig = Figure(size = (1000, 500))
ax1 = Axis(fig[1:2,1])
ax2 = Axis(fig[1,2], title = "Sample rate = $fps hz", ylabel = "Fluorescence ($(data2P.chUnits[1]))")
ax3 = Axis(fig[2,2], title = "Sample rate = $(1/dataIC.dt) hz", ylabel = "Membrane Voltage ($(dataIC.chUnits[1]))")
hm1 = heatmap!(ax1, xlims, ylims, zproj, colormap = Reverse(:speed), colorrange = (0.001, 0.01))
lines!(ax2, t, fluo, color = :green)

experimentplot!(ax3, dataIC, channel = 1)
ticker1 = vlines!(ax2, [0.0], color = :black)
ticker2 = vlines!(ax3, [0.0], color = :black)
Colorbar(fig[3, 1], hm1, vertical = false, label = "Fluorescence (px)")

# Animate the data ___________________________________________________________________#
record(fig, save_fn, enumerate(t), framerate = fps*10) do (i, t) 
    println(t)
    hm1[3] = mov[:,:,i]
    ticker1[1] = [t]
    ticker2[1] = [t]
end
display(fig)



# Open and prepare the data______________________________________________________________#
fps = 1.45
data2P = readImage(filepath; sampling_rate = fps, chName = "525 nm");
px_x, px_y = data2P.HeaderDict["framesize"]
xlims = data2P.HeaderDict["xrng"]; #Extract the x domain
ylims = data2P.HeaderDict["yrng"]; #Extract the y domain
time_2P = data2P.t; #Extract the t domain
mov = get_all_frames(data2P) #get all frames as a 3D array
fluo = mean(mov, dims = (1,2))[1,1,:]
zproj = maximum(mov, dims = 3)[:,:,1] #Can we save this directly? 

#%% Point to the file location ____________________________________________________________#
drive = "D:/"
root = "Data/Patching/"
folder = "2024_03_11_VglutGC6/Cell3/"
file = "24311015.abf"
filepath = joinpath(drive, root, folder, file)
# Open and prepare the data______________________________________________________________#
dataIC = readABF(filepath); 
downsample!(dataIC, 1000.0); #Go through these functions much more extensively
time_IC = dataIC.t

#%% Create the plot_______________________________________________________________________#
fig = Figure(size = (1000, 500))
ax1 = Axis(fig[1:2,1])
ax2 = Axis(fig[1,2], title = "Sample rate = $fps hz", ylabel = "Fluorescence ($(data2P.chUnits[1]))")
ax3 = Axis(fig[2,2], title = "Sample rate = $(1/dataIC.dt) hz", ylabel = "Membrane Voltage ($(dataIC.chUnits[1]))")
hm1 = heatmap!(ax1, xlims, ylims, zproj, colormap = Reverse(:speed), colorrange = (0.001, 0.01))
lines!(ax2, t, fluo, color = :green)
lines!(ax3, dataIC.t, dataIC.data_array[1,:, 1], color = :black)

ticker1 = vlines!(ax2, [0.0], color = :black)
ticker2 = vlines!(ax3, [0.0], color = :black)
Colorbar(fig[3, 1], hm1, vertical = false, label = "Fluorescence (px)")

# Animate the data ___________________________________________________________________#
record(fig, "D:/Data/Analysis/2024_03_11_VGlutGC6_Cell3_2PIC.mp4", 2:size(data,2), framerate = fps*10) do i 
    hm1[3] = mov[:,:,i]
    ticker1[1] = [t[i]]
    ticker2[1] = [t[i]]
end
display(fig)