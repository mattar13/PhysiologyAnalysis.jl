#%% Using the Julia interface
using Revise, ePhys 
using Pluto
run_experiment_analysis()
#run_filter_determination()

#%% Try to save ABF
using Revise
using ePhys
import ePhys.saveABF
import ePhys.Experiment
using PyPlot
Revise.track(ePhys, "src/Readers/ABFReader/ABFReader.jl")
file_open = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Cones\2019_07_23_WT_P14_m1\Cones\Drugs\Green\nd0.5_1p_1ms\19723190.abf"
file_save = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Cones\2019_07_23_WT_P14_m1\Cones\Drugs\Green\nd0.5_1p_1ms\test.abf"
data = readABF(file_open, channels=-1) #This is necessary for saving
saveABF(data, file_save)

#%% Open Pauls files
ePhys.__init__()
using DataFrames, Query, XLSX
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)

#%% Can we open CSV files? 
using Revise, ePhys
using PyPlot
import ePhys.plot_experiment
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Human\Donor07_RPeriphery_1.csv"

data = readCSV(file)
#%%
fig,ax = plt.subplots(1)
plot_experiment(ax, data)
ax.set_xlabel("Time (s)")
ax.set_xlim(0.0, 3.0)
ax.set_ylabel("Response (μV)")
ax.set_ylim(-350, 200)
#%%
a = -minimum(data, dims = 2)
b = maximum(data, dims = 2)
ax[2].scatter(a, b.-a)

#%% Test out some git stuff here
println("This is the master branch")