using ElectroPhysiology, PhysiologyAnalysis

#Opening and reanalyzing a previously analyzed file
data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2023_05_24 Critical Timepoints in retinoschisis - Unknown Journal" #The data root

#%% Opening a previously existant file
data_file = joinpath(data_root, "data_analysis.xlsx") #This is the main file we can use. The root may change
dataset = openDataset(data_file)

#%% opening a file from new
files = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_12_20_WTAdult\Mouse1_WT_Adult"|>parseABF
dataset = createDataset(files, verbose = true)

dataset = runTraceAnalysis(dataset, verbose = true)
dataset = runExperimentAnalysis(dataset)
dataset = runConditionsAnalysis(dataset)
backupDataset(data_file)
saveDataset(dataset, data_file)


dataset = createDataset(files, verbose = true)
