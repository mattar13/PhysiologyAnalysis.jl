using Revise
using ElectroPhysiology
using PhysiologyAnalysis



#%%
using PyPlot
PyPlot.pygui(true)
import PhysiologyAnalysis: readABF, parseABF
import PhysiologyAnalysis: DataFrame
using DataFrames, Query, XLSX

data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis" #The data root |> parseABF
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\Retinoschisis\\data_analysis.xlsx"
dataset = openDataset(datafile)
dataset["CONDITIONS"] = summarize_data(dataset)
dataset["STATS"] = dataset_statistics(dataset)
backupDataset(datafile)
saveDataset(dataset, datafile)

plot_dataset_fits(dataset, normalize = false, condition = "NoDrugs")
plot_dataset_vals(dataset)

#test quick plot
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2023_02_23_MattR141C\Mouse2_Adult_R141C\BaCl_LAP4\Rods"
data = readABF(file |> parseABF) |> data_filter
plot_experiment(data, xlims = (-0.25, 1.0))

#%% Make a new function for adding photon numbers to a flash stimulus
using PhysiologyAnalysis
file_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_04_21_a13MelCreAdult\Mouse2_Adult_WT"
files_A = joinpath(file_root, "BaCl_LAP4\\Rods") |> parseABF