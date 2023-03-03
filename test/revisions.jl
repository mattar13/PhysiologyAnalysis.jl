#%% Section 1. Revision of interface
using Revise
using PhysiologyAnalysis
import PhysiologyAnalysis: readABF, parseABF
import PhysiologyAnalysis: DataFrame
#%% Add the ability to flag files and remove them from the analysis

testRoot = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2023_02_23_MattR141C\Mouse1_Adult_R141C"
dataset = createDataset(testRoot)
dataset["STATS"] = dataset_statistics(dataset)
dataset["EXPERIMENTS"]