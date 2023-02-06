using Revise, ePhys
using DataFrames, Query, XLSX
import ePhys: run_A_wave_analysis, run_B_wave_analysis, run_G_wave_analysis
exp_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_07_11_Adult\Mouse3_Adult_RS1KO"
experiment_paths = exp_root |> parseABF

all_files = createDatasheet(experiment_paths; filename = savefile_name)

res_A = run_A_wave_analysis(all_files)