#=[Import packages]============================================================#
using ElectroPhysiology, PhysiologyAnalysis
using Statistics

#=[Point to filenames]=========================================================#
data_fn = "D:/Data/Patching/2024_03_11_VglutGC6/Cell2/24311005.abf"
save_fn = "PatchPlot.png"

#=[Open the data]==============================================================#
data = readABF(data_fn); 
downsample!(data, 1000.0); #Go through these functions much more extensively
t = data.t

#=[Do some analysis]===========================================================#
threshold = calculate_threshold(data, Z = 2)
timestamps = get_timestamps(data, Z = 2)
#durations, intervals = extract_interval(timestamps[1,1])
# Eventually need to figure this out, but not today
# max_interval_algorithim(timestamps[1,1], SPBmin = 1, verbose = true)

#=[Plot data]==================================================================#
fig, axs = experimentplot(data)
save(save_fn, fig)