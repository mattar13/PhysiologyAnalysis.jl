using Revise
using ElectroPhysiology 
using PhysiologyAnalysis
using Pkg; Pkg.activate("test")
using GLMakie
using FileIO, ImageView, Images, ImageFiltering
using Statistics
import GLMakie.Axis

#%% Open the data Set the domain and parameters
domain_x = (xmin, xmax) = (0.0, 0.3445)
domain_y = (ymin, ymax) = (0.0, 0.3445)

data2P_fn = raw"G:\Data\Calcium Imaging\2024_05_23_MORF_ChATCre\ca_img_5011.tif"
data2P = readImage(data2P_fn)
deinterleave!(data2P) #This seperates the movies into two seperate movies
t = data2P.t
img_arr = get_all_frames(data2P)
ca_img = img_arr[:,:,:,1]
ca_trace = project(data2P, dims = (1,2))[1,1,:,1]
maximum(ca_trace)
minimum(ca_trace)

green_ch = project(data2P, dims = 3)[:,:,1,1]
red_ch = project(data2P, dims = 3)[:,:,1,2] #Take the mean of the red channel to get the cell marker
red_ch[findall(red_ch .< 0.005)] .= NaN #Set the mask

xlims = LinRange(ymin, xmax, size(img_arr,1))
ylims = LinRange(ymin, ymax, size(img_arr,2))
dx = xlims[2]-xlims[1]
dy = ylims[2]-ylims[1]

# Calculate the delta F over F
voxel_z = 50
background = roll_mean(ca_img; voxel_z = voxel_z, mode = "reflect")
background_trace = mean(background, dims = (1,2))[1,1,:]

delta_f_f = deltaF_F(ca_img; voxel_z = voxel_z, mode = "reflect")
delta_f_f_mean = maximum(delta_f_f, dims = 3)[:,:,1]
delta_f_f_trace = mean(delta_f_f, dims = (1,2))[1,1,:]

#%% Plot the figure
fig1 = Figure(size = (1000, 500))
ax1 = GLMakie.Axis(fig1[1:2,1], title = "Dual Experiment", aspect = 1.0)
#heatmap!(ax1, xlims, ylims, green_ch, colormap = Reverse(:speed), colorrange = (0.0, maximum(green_ch)/10))
#heatmap!(ax1, xlims, ylims, delta_f_f_mean, colormap = Reverse(:speed), colorrange = (0.0, maximum(delta_f_f_mean)/10))

hm1 = heatmap!(ax1, xlims, ylims, delta_f_f[:,:,1], colormap = Reverse(:speed), colorrange = (0.0, maximum(delta_f_f)/10))#, alpha = 0.25)
heatmap!(ax1, xlims, ylims, red_ch, colormap = Reverse(:amp), colorrange = (0.0, maximum(red_ch)), alpha = 0.4)

ax2a = GLMakie.Axis(fig1[1,2])
lines!(ax2a, t, ca_trace)
lines!(ax2a, t, background_trace, color = :red)
ticker1 = vlines!(ax2a, [0.0], color = :black)

ax2b = GLMakie.Axis(fig1[2,2])
lines!(ax2b, t, delta_f_f_trace)
ticker2 = vlines!(ax2b, [0.0], color = :black)

display(fig1)

#%% Save the figure once you are ready
save_fn = "test/delta_f_f.png"
save(save_fn, fig1)

#%% Run the animation
save_fn = "test/delta_f_f.mp4"
fps = 1/data2P.dt
GLMakie.record(fig1, save_fn, enumerate(data2P.t), framerate = 5fps) do (i, t) 
     println(t)
     hm1[3] = delta_f_f[:,:,i]
     ticker1[1] = [t]
     ticker2[1] = [t]
 end