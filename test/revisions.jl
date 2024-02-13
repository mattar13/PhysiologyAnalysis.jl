using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Pkg; Pkg.activate("test")
PhysiologyAnalysis.__init__()
using FileIO, ImageView, Images
using GLMakie
import GLMakie.Axis
using Flux

#%% 
root = raw"F:\Data\Calcium Images\2024_01_31_FRMD7_Cal590"
file = "2024_01_31_ROI004.tif"
fn = joinpath(root, file)
data = readImage(fn)
mov = get_all_frames(data)

xsize, ysize, zsize = size(mov)
xlims = 1:xsize
ylims = 1:ysize
zlims = 1:zsize
#Create an average for the frames
zavg = (sum(mov, dims = 3) ./ size(mov, 3))[:,:,1]

#Minimize the zavg
threshold = 0.5
zproj = (sum(data.data_array, dims = 1) ./ size(data, 1))[1, :, 1]

fig = Figure()
ax1 = GLMakie.Axis(fig[1,1])
heatmap!(ax1, xlims, ylims, zavg, colormap = Reverse(:algae))
display(fig)
save("test.png", fig)
save("test_img.png", zavg)
#%%
recordROI(data, [1,2,3,4], 1)
getROImask(data, 1)

ROI1 = getROIarr(data, 1)

#%%
fig = GLMakie.Figure(size = (1500, 400))
ax11 = GLMakie.Axis(fig[1,1])
ax12 = GLMakie.Axis3(fig[1,2])
ax21 = GLMakie.Axis(fig[2,1])
ax22 = GLMakie.Axis3(fig[2,2])
ax3 = GLMakie.Axis(fig[3,1])

heatmap!(ax11, xlims, ylims, zavg, colormap = Reverse(:algae), colorrange = (0.0, 0.05))
surface!(ax12, xlims, ylims, zavg, colormap = Reverse(:algae), colorrange = (0.0, 0.05))
hm = heatmap!(ax21, xlims, ylims, mov[:,:,1], colormap = Reverse(:algae), colorrange = (0.0, 0.05))
sf = surface!(ax22, xlims, ylims, mov[:,:,1], colormap = Reverse(:algae), colorrange = (0.0, 0.05))
lines!(ax3, zproj)
ln_spot = lines!(ax3, fill(0.0, 2), [minimum(zproj), maximum(zproj)], color = :black)

record(fig, "ANIMATION.mp4", 2:size(mov,3)) do i 
    println(i)
    frame = mov[:,:,i]
    ln_spot[1] = fill(i, 2)
    hm[3] = frame
    sf[3] = frame
end
display(fig)