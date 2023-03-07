#%% Section 1. Revision of interface
using Revise
using PhysiologyAnalysis
using PyPlot
PyPlot.pygui(true)
import PhysiologyAnalysis: readABF, parseABF
import PhysiologyAnalysis: DataFrame
using DataFrames, Query, XLSX

#enter the topmost file root here (date, mouse info)
exp_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2021_09_12_RS1KO-13\Mouse1_P13_RS1KO"
channels = ["Vm_prime"]; #Specify which channels you are using
photoreceptors = "Rods";
dataset = createDataset(exp_root; t_post = 2.0, verbose = false);
trace_A = matchExperiment(dataset["TRACES"], (Condition = "BaCl", Photoreceptor = "Rods"))
exp_A = matchExperiment(dataset["EXPERIMENTS"], (Condition = "BaCl", Photoreceptor = "Rods"))
flagExperiment!(dataset["EXPERIMENTS"], exp_A)

#%% Make a plot data summary function
plot_data_summary(trace_A, exp_A)

#%% Add the ability to flag files and remove them from the analysis
data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis" #The data root |> parseABF
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\Retinoschisis\\data_analysis.xlsx"
all_files = openDataset(datafile, sheetnames = "ALL_FILES")
dataset = runTraceAnalysis(all_files)
dataset["STATS"] = dataset_statistics(dataset)
backupDataset(datafile)
saveDataset(dataset, datafile)