#%% Load Packages _________________________________________________________________________#
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyAnalysis
Pkg.activate("test")
using GLMakie
import GLMakie.Axis
using FileIO, ImageView, Images
using Statistics

#%% Point to the file location ____________________________________________________________#
drive = "D:/"
root = "Data/Calcium Imaging/"
folder = "2024_03_11_VglutGC6/Cell2/"
file = "freeRec004.tif"
filepath = joinpath(drive, root, folder, file);

# Open and prepare the data______________________________________________________________#
fps = 1.45
data = readImage(filepath; sampling_rate = fps);
px_x, px_y = data.HeaderDict["framesize"]
xlims = data.HeaderDict["xrng"]; #Extract the x domain
ylims = data.HeaderDict["yrng"]; #Extract the y domain
t = data.t; #Extract the t domain

mov = get_all_frames(data); #get all frames as a 3D array
fluo = mean(mov, dims = (1,2))[1,1,:]
zproj = maximum(mov, dims = 3)[:,:,1] #Can we save this directly? 

# Plot the data _________________________________________________________________________#
fig = Figure(size = (500, 750))
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
rowsize!(fig.layout, 2, Relative(1/4))
hm1 = heatmap!(ax1, xlims, ylims, zproj, colormap = Reverse(:speed), colorrange = (0.0, 0.01))
lines!(ax2, t, fluo)
ticker = vlines!(ax2, [0.0])
display(fig)

#%% Animate the Imaging video ___________________________________________________________#
record(fig, "D:/Data/Analysis/ca_img.mp4", 2:size(data,2), framerate = fps*10) do i 
    println(i)
    hm1[3] = mov[:,:,i]
    ticker[1] = [t[i]]
end