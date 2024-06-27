using Revise
using Dates
using ElectroPhysiology 
using PhysiologyAnalysis
using Pkg; Pkg.activate("test")

# Can we do an analysis for the IV curves
data_fn = raw"D:\Data\Patching\2024_06_22_VCGC6_P8\Cell4\24622021.abf"

data = readABF(data_fn); #Open the data
V_HOLD = extract_timepoint(data, channel = 2)
I_CM = calculate_peak(data)
I_RIN = extract_timepoint(data)
Rs, Rin = calculate_resistance(data)
Cm = calculate_capacitance(data)

PhysiologyAnalysis.__init__()

#%% When we get to the dataframe creation we can fill in details here
using XLSX, DataFrames, Query
using FileIO, Images, ImageView, ImageMagick
import .PhysiologyAnalysis.traverse_root

using Statistics
using PhysiologyPlotting, GLMakie
import PhysiologyPlotting.experimentplot
using CurveFit
import CurveFit.curve_fit

#=[Point to filenames]=========================================================#
data_fn = raw"D:\Data\Patching\2024_06_22_VCGC6_P8\Cell4\24622021.abf"

#=[Open the data]==============================================================#
data = readABF(data_fn); 
#downsample!(data, 100.0); #Go through these functions much more extensively

#Extract the IV
selected_idx = findfirst(data.t .>= 0.5)
IV_vals = data[:, selected_idx, :]
selected_idx = findfirst(data.t .>= 0.1317)

#Calculate the series resistance
VC_hold = data[:, selected_idx, 2]
IC_hold = data[:, selected_idx, 1]
IV_fit = CurveFit.curve_fit(LinearFit, IC_hold, VC_hold)
series_resistance = IV_fit.coefs[1]

Ic_fit = LinRange(-1000, 3000, 1000)
Vc_fit = IV_fit.(Ic_fit)

#%%=[Plot data]==================================================================#
fig = Figure(size = (900, 400))
ax1a = Axis(fig[1,1], ylabel = "Current (pA)")
ax2a = Axis(fig[2,1], xlabel = "Time (s)", ylabel = "Voltage (mV)")
linkxaxes!(ax1a, ax2a)
ax1b = Axis(fig[1:2,2], xlabel = "Voltage (mV)", ylabel = "Current (pA)")
experimentplot!(ax1a, data, channel = 1)
experimentplot!(ax2a, data, channel = 2)
scatterlines!(ax1b, IV_vals[:, 2], IV_vals[:, 1], markersize = 10.0, color = :black)
#scatter!(ax1b, VC_hold, IC_hold, markersize = 10.0, color = :red)
#lines!(ax1b, Vc_fit, Ic_fit, color = :red)
hlines!(ax1b, [0.0], color = :black)
vlines!(ax1b, [0.0], color = :black)
display(fig)

#%% Open all dates for the imaging files
img_dir = raw"D:\Data\Calcium Imaging"
patch_dir = raw"D:\Data\Patching"
filename = raw"D:\Data\Analysis\data_analysis.xlsx"

#Don't often need to do this (eventually we need to add an update button)
#dataset = create2PDataSheet(img_dir, patch_dir; verbose = true)

PhysiologyAnalysis.__init__()
dataset = open2PDataSheet(filename)
all_datasheet = dataset["All Files"]
img_datasheet = dataset["TIF Files"]
patch_datasheet = dataset["ABF Files"] 
dataset["Paired Files"] = pair_experiments(patch_datasheet, img_datasheet)
save2PDataSheet(filename, dataset)
