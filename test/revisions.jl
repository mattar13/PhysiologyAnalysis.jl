#%% Section 1. Revision of interface
using Revise
using PhysiologyAnalysis
import PhysiologyAnalysis: readABF, parseABF
import PhysiologyAnalysis: DataFrame
#%% Add the ability to flag files and remove them from the analysis
datafile = raw"C:\Users\mtarc\OneDrive - The University of Akron\Projects\Retinoschisis\data_analysis.xlsx" #The data root
dataset = openDataset(datafile, sheetName="all", typeConvert = true)
dataset["STATS"] = dataset_statistics(dataset["EXPERIMENTS"])
saveDataset(dataset, datafile)
backupDataset(datafile)