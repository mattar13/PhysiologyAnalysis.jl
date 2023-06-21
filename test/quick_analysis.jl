using ElectroPhysiology
using DataFrames, Query, XLSX
using PhysiologyAnalysis

#Opening and reanalyzing a previously analyzed file
data_file = save_loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2023_05_24 Critical Timepoints in retinoschisis - Unknown Journal\data_analysis.xlsx" #The data root
dataset = openDataset(data_file)
backupDataset(data_file) #Backup the datafile

#%% opening a file from new
files = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2019_03_12_AdultWT"|>parseABF
save_loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2023_05_24 Critical Timepoints in retinoschisis - Unknown Journal\2019_03_12_Analysis.xlsx"
dataset = createDataset(files, verbose = true)

#%% Conducting the analysis
dataset = runTraceAnalysis(dataset, verbose = true)
dataset = runExperimentAnalysis(dataset, verbose = true)
dataset = runConditionsAnalysis(dataset)
dataset = runStatsAnalysis(dataset)

#%% Saving the data
saveDataset(dataset, save_loc)