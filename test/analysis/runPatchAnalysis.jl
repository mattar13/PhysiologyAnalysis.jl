#%% Load Packages _________________________________________________________________________#
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyAnalysis
Pkg.activate("test")
using GLMakie
using Statistics

#%% Point to the file location ____________________________________________________________#
drive = "D:/"
root = "Data/Patching/"
folder = "2024_03_11_VglutGC6/Cell3/"
file = "24311012.abf"
filepath = joinpath(drive, root, folder, file)

# Open and prepare the data______________________________________________________________#
data = readABF(filepath); 
downsample!(data, 1000.0); #Go through these functions much more extensively
t = data.t

# Do some analysis _____________________________________________________________________#
threshold = calculate_threshold(data, Z = 2)
timestamps = get_timestamps(data, Z = 2)
#durations, intervals = extract_interval(timestamps[1,1])
# Eventually need to figure this out, but not today
# max_interval_algorithim(timestamps[1,1], SPBmin = 1, verbose = true)

#%% Plot the data _________________________________________________________________________#
fig = Figure(title = "Test")
ax1 = Axis(fig[1,1], title = filepath,
     ylabel = "$(data.chNames[1]) ($(data.chUnits[1]))")
ax2 = Axis(fig[2,1], xlabel = "Time (ms)", ylabel = "$(data.chNames[2]) ($(data.chUnits[2]))")

for i in 1:size(data,1)
     lines!(ax1, data.t, data.data_array[i,:, 1], colormap = :berlin, color = [i], colorrange = (1, size(data,1)))
     lines!(ax2, data.t, data.data_array[i,:, 2], colormap = :berlin, color = [i], colorrange = (1, size(data,1)))
end

display(fig)
save("D:/Data/Analysis/test.png", fig)