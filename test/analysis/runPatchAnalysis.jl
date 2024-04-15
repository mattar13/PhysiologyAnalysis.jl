#=[Import packages]============================================================#
using Pkg; Pkg.activate(".") 
using ElectroPhysiology, PhysiologyAnalysis
using Statistics
Pkg.activate("test")
using PhysiologyPlotting
using GLMakie
import PhysiologyPlotting.experimentplot
#=[Point to filenames]=========================================================#
data_fn = "D:/Data/Patching/2024_02_12_WT/Cell1/24212004.abf"
save_fn = "D:/Data/Analysis/2024_02_12_WT_Cell1_IC.png"

#=[Open the data]==============================================================#
data = readABF(data_fn); 
downsample!(data, 1000.0); #Go through these functions much more extensively
data.HeaderDict
#=[Plot data]==================================================================#
fig, axs = experimentplot(data)
save(save_fn, fig)

#=[Do some analysis]===========================================================#
threshold = calculate_threshold(data, Z = 2)
timestamps = get_timestamps(data, Z = 2)
#durations, intervals = extract_interval(timestamps[1,1])
# Eventually need to figure this out, but not today
# max_interval_algorithim(timestamps[1,1], SPBmin = 1, verbose = true)
