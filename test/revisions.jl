#%% Section 1. Revision of interface
using Revise
using ElectroPhysiology
using PhysiologyAnalysis

import PhysiologyAnalysis: readABF, parseABF
#%% Add the ability to flag files and remove them from the analysis
data_file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Projects\Retinoschisis\data_analysis.xlsx" #The data root
dataset = openDataset(data_file, sheetName="all", typeConvert = true)