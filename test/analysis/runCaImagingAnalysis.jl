#%% Load Packages _________________________________________________________________________#
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyAnalysis
Pkg.activate("test")
using PhysiologyPlotting, GLMakie
import GLMakie.Axis
using FileIO, ImageView, Images
using Statistics

#=[Point to filenames]=========================================================#
data_fn = "D:/Data/Calcium Imaging/2024_02_12_WT/Cell1/cell1001.tif"
save_fn = "D:/Data/Analysis/2024_02_12_WTCell1_2P.mp4"

#[Open and prepare the data]______________________________________________________________#
fps = 1.45
data = readImage(data_fn);
px_x, px_y = data.HeaderDict["framesize"]
xlims = data.HeaderDict["xrng"]; #Extract the x domain
ylims = data.HeaderDict["yrng"]; #Extract the y domain
t = data.t #Extract the t domain
#Every other frame is a different channel
t = t[1:2:end]

mov = get_all_frames(data)[:,:,1:2:end]; #get all frames as a 3D array
fluo = mean(mov, dims = (1,2))[1,1,:]
zproj = maximum(mov, dims = 3)[:,:,1] #Can we save this directly? 

# Plot the data _________________________________________________________________________#
fig = Figure(size = (500, 600))
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
rowsize!(fig.layout, 2, Relative(1/4))
hm1 = heatmap!(ax1, xlims, ylims, zproj, colormap = Reverse(:speed), colorrange = (0.0, 0.01))
lines!(ax2, t, fluo)
ticker = vlines!(ax2, [0.0])

# Animate the Imaging video ___________________________________________________________#
record(fig, save_fn, 2:length(t), framerate = fps*10) do i 
    println(i)
    hm1[3] = mov[:,:,i]
    ticker[1] = [t[i]]
end