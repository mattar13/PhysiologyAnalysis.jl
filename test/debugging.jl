using Revise
using ePhys
using PyPlot
using DataFrames, Query, XLSX
ePhys.__init__()

#%% Why is the data from wildtype looking so weird
file_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Paul\2019_09_24_WT-30\Mouse1_Adult_WT\BaCl_LAP4\Rods"
files = file_root |> parseABF
df = createDatasheet(files; filename = nothing)
#%% Open the datafile
data = data_filter(readABF(files, channels = ["Vm_prime", "IN 7"]), 
     scale = 1.0, #Keep the units in mV
     filter_channels = "Vm_prime", freq_stop = 40.0, t_pre = 0.25, t_post = 3.75, remove_global_drift = false
)
#How can we downsample the data correctly

downsample!(data, 1/0.001)

fig, ax = plt.subplots(2)
plot_experiment(ax[1], data, channels = 1)
plot_experiment(ax[2], data, channels = 2)

#%% Open Pauls files
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)