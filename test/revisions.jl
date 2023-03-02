#%% Section 1. Revision of interface
using Revise
using ElectroPhysiology
using PhysiologyAnalysis

import PhysiologyAnalysis: readABF, parseABF
#%% Add the ability to flag files and remove them from the analysis
exp_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2023_02_23_MattR141C\Mouse1_Adult_R141C\BaCl_LAP4\Rods"
dataset = createDataset(exp_root)