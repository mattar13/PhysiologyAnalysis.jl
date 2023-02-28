#using Revise
using PhysiologyAnalysis
using DataFrames, Query, XLSX
import PhysiologyAnalysis.restructure_filesystem
#%% New file location
#prev_paths = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Rods" |> parseABF
#prev_paths = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\NotDetermined" |> parseABF
#prev_paths = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Cones" |> parseABF
prev_paths1 = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\JGP_Files\Cones_photopic" |> parseABF
prev_paths2 = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\JGP_Files\Rods_scotopic" |> parseABF
prev_paths = [prev_paths1..., prev_paths2...]
rec_path = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\JGPRestructured"
restructure_filesystem(prev_paths, rec_path)
rec_path |> parseABF