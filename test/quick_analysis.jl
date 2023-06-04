using ElectroPhysiology, PhysiologyAnalysis

data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis" #The data root
all_files = data_root |> parseABF
dataset["ALL_FILES"]
dataset = createDataset(all_files, verbose = true)
dataset = runTraceAnalysis(dataset, verbose = true)
dataset = runExperimentAnalysis(dataset, verbose = true)
datafile = raw"C:\Users\mtarc\OneDrive - The University of Akron\Projects\Retinoschisis\data_analysis.xlsx"

#%%
using DataFrames, Query
using ElectroPhysiology, PhysiologyAnalysis
datafile = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2023_05_24 Critical Timepoints in retinoschisis - Unknown Journal\data_analysis.xlsx" #The data root
dataset = openDataset(datafile)
size(dataset["EXPERIMENTS"])
dataset["EXPERIMENTS"] |> @filter(_.Month == "12" && _.Date == "20" && _.Number == "1") |> DataFrame
dataset["EXPERIMENTS"] |> @filter(_.Age == "Adult" && _.Genotype == "WT" && _.Month == "04") |> DataFrame

dataset = runExperimentAnalysis(dataset)
saveDataset(dataset, datafile)