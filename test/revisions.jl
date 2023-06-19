using Dates
using Revise
using PhysiologyAnalysis
using ElectroPhysiology
using DataFrames, Query

#%% Make a new function for adding photon numbers to a flash stimulus
using PhysiologyAnalysis
file_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_04_21_a13MelCreAdult\Mouse2_Adult_WT"
files_A = joinpath(file_root, "BaCl_LAP4\\Rods") |> parseABF