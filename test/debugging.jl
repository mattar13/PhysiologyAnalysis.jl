#%% Set up plotting so each trace can be colored based on a gradient
using Revise
using ePhys
using PyPlot
import ePhys.plot_experiment
paths = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Data\\ERG\\Gnat\\2020_12_12_gnat14\\Mouse1_P14_GNAT-KO\\BaCl_LAP4\\525Green" |> parseABF
data = readABF(paths, channels = ["Vm_prime"]) |> data_filter

#%%
fig,ax = plt.subplots(1)
plot_experiment(data, color = :black)
ePhys.is_cmap(:jet)
cmap = plt.get_cmap(:jet)
cmap(0.5)

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

using Revise

Revise.add_file