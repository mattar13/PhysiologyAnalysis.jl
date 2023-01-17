using ePhys
using DataFrames, Query, XLSX

#%% Run the data analysis for the retinoschisis dataset
data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis" #The data root
datafile = raw"C:\Users\mtarc\OneDrive - The University of Akron\Projects\Retinoschisis\data_analysis.xlsx"
all_files = data_root |> parseABF
dataset = openDatasheet(datafile, sheetName="all")
updateDatasheet(datafile, all_files)
runAnalysis(datafile)

#%% Run the data analysis for the retinoschisis dataset
data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Gnat" #The data root
datafile = raw"C:\Users\mtarc\OneDrive - The University of Akron\Projects\GNAT\data_analysis.xlsx"
all_files = data_root |> parseABF
createDatasheet(all_files; filename = datafile, verbose = true)
dataset = openDatasheet(datafile, sheetName="all")
updateDatasheet(datafile, all_files)
runAnalysis(datafile)
