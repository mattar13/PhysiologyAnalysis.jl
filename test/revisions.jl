using Revise
using ElectroPhysiology
import ElectroPhysiology as EP
using PhysiologyAnalysis
ElectroPhysiology.__init__()
using Pkg; Pkg.activate("test")
using FileIO, ImageView, Images
ElectroPhysiology.Images
ElectroPhysiology.ImageView

using GLMakie
using Statistics

root = raw"F:\Data\Calcium Images\2024_01_31_FRMD7_Cal590"
file = "2024_01_31_ROI004.tif"
mask_file = "2024_01_31_ROI004_MASK.png"

fn = joinpath(root, file)
mask_fn = joinpath(root, file)
loadROIfn!(mask_fn, data)
data = readImage(fn)
mov = get_all_frames(data) #get all frames as a 3D array
zproj = mean(mov, dims = 3)[:,:,1] #Can we save this directly? 

#%% Figure out how to load a mask image as a ROI array
ROI_mask = getROImask(data)

px_x, px_y = data.HeaderDict["framesize"]
xlims = 1:px_x
ylims = 1:px_y
zlims = 1:size(data,2)
#%% Lets make a plot that shows a few of the ROIs
fig = GLMakie.Figure(size = (1000, 500))
#Plot the projected image and the mask
ax11 = GLMakie.Axis(fig[1,1])
ax12 = GLMakie.Axis(fig[1,2])
heatmap!(ax11, xlims, ylims, zproj, colormap = Reverse(:algae))
heatmap!(ax12, xlims, ylims, ROI_mask, colormap = Reverse(:algae))

# Extract the ROI movie array
ROI1_mov = getROIarr(data, 21)
ROIall_mov = getROIarr(data)
ax21 = GLMakie.Axis(fig[2,1])
ax22 = GLMakie.Axis(fig[2,2])
hm1 = heatmap!(ax21, xlims, ylims, mov[:,:,1], colormap = Reverse(:algae), colorrange = (0.0, 0.05))
hm2 = heatmap!(ax22, xlims, ylims, ROIall_mov[:,:,1], colormap = Reverse(:algae), colorrange = (0.0, 0.05))

#Convert the centroid of each cell coordinate to an X,Y coordinate
display(fig)
#%%
record(fig, "ANIMATION.mp4", 2:size(data,2)) do i 
     println(i)
     hm1[3] = mov[:,:,i]
     hm2[3] = ROIall_mov[:,:,i]
end
display(fig)