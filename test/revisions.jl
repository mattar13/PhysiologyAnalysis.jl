using Revise
using Dates
using ElectroPhysiology 
using PhysiologyAnalysis
using Pkg; Pkg.activate("test")

#=[Point to filenames]=========================================================#
data_fn = raw"G:\Data\Patching\2024_06_22_VCGC6_P8\Cell4\24622021.abf"

#=[Open the data]==============================================================#
PhysiologyAnalysis.__init__()

data = readABF(data_fn); 
V_HOLD = extract_timepoint(data, channel = 2)
I_CM = calculate_peak(data)
I_RIN = extract_timepoint(data)
Rs, Rin, V_M, V_HOLD = calculate_resistance(data)
Cm = calculate_capacitance(data)
τ = Rs*Cm
Quality = Rs/Rin * 100

#%% =[Plot data]==================================================================#
using GLMakie, PhysiologyPlotting
fig = Figure(size = (900, 400))
ax1a = GLMakie.Axis(fig[1,1], ylabel = "Current (pA)")
ax2a = GLMakie.Axis(fig[2,1], xlabel = "Time (s)", ylabel = "Voltage (mV)")
linkxaxes!(ax1a, ax2a)
ax1b = GLMakie.Axis(fig[1:2,2],  xlabel = "Voltage (mV)", ylabel = "Current (pA)")
experimentplot!(ax1a, data, channel = 1, xlims = (0.0, 1.0))
experimentplot!(ax2a, data, channel = 2, xlims = (0.0, 1.0))
#xlims!(ax1a, (0.1312, 0.6312))
#xlims!(ax2a, (0.1312, 0.6312))

scatterlines!(ax1b, V_HOLD, I_RIN, markersize = 10.0, color = :black)
scatterlines!(ax1b, V_HOLD, I_CM, markersize = 10.0, color = :red)

hlines!(ax1b, [0.0], color = :black)
vlines!(ax1b, [0.0], color = :black)
display(fig)

#%% When we get to the dataframe creation we can fill in details here
using XLSX, DataFrames, Query
using FileIO, Images, ImageView, ImageMagick
import .PhysiologyAnalysis.traverse_root

#%% Open all dates for the imaging files
img_dir = raw"D:\Data\Calcium Imaging"
patch_dir = raw"D:\Data\Patching"
filename = raw"D:\Data\Analysis\data_analysis.xlsx"

#Don't often need to do this (eventually we need to add an update button)
#dataset = create2PDataSheet(img_dir, patch_dir; verbose = true)

dataset = open2PDataSheet(filename)
all_datasheet = dataset["All Files"]
img_datasheet = dataset["TIF Files"]
patch_datasheet = dataset["ABF Files"] 

IV_analysis!(dataset)
size(patch_datasheet)
dataset = pair_experiments!(dataset)
save2PDataSheet(filename, dataset)
