#%% Load Packages _________________________________________________________________________#
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyAnalysis
Pkg.activate("test")
using GLMakie

#%% Point to the file location ____________________________________________________________#
root = "F:/Data/Patching/"
folder = "2024_03_11_VglutGC6/Cell3"
file = "24311012.abf"
filepath = joinpath(root, folder, file)

# Open and prepare the data______________________________________________________________#
data = readABF(filepath); 
downsample!(data, 1000.0); #Go through these functions much more extensively
t = data.t

#%% Plot the data _________________________________________________________________________#
fig = Figure()
ax1 = Axis(fig[1,1], ylabel = "$(data.chNames[1]) ($(data.chUnits[1]))")
ax2 = Axis(fig[2,1], xlabel = "Time (ms)", ylabel = "$(data.chNames[2]) ($(data.chUnits[2]))")

for i in 1:size(data,1)
     lines!(ax1, data.t, data.data_array[i,:, 1], colormap = :berlin, color = [i], colorrange = (1, size(data,1)))
     lines!(ax2, data.t, data.data_array[i,:, 2], colormap = :berlin, color = [i], colorrange = (1, size(data,1)))
end
display(fig)

#%% Do some analysis _____________________________________________________________________#