#%% Section 1. Revision of interface
using Revise
using PhysiologyAnalysis
import PhysiologyAnalysis: readABF, parseABF
import PhysiologyAnalysis: DataFrame

#%% Add the ability to flag files and remove them from the analysis

data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis" #The data root |> parseABF
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\Retinoschisis\\data_analysis.xlsx"
all_files = openDataset(datafile, sheetnames = "ALL_FILES")
dataset = runTraceAnalysis(all_files)
dataset["STATS"] = dataset_statistics(dataset)
saveDataset(dataset, all_files)