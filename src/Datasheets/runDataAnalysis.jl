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
createDatasheet(all_files; filename=datafile, verbose=true)
dataset = openDatasheet(datafile, sheetName="all")
updateDatasheet(datafile, all_files)
runAnalysis(datafile)

#%% Run the data analysis for the JGP files
using Dates
Dates.now()
yyyy = year(Dates.now())
mm = month(Dates.now())
dd = day(Dates.now())
data_root  = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\JGP_Files\Traces"
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\GNAT\\$(yyyy)_$(mm)_$(dd)_cone_data_analysis.xlsx"
all_files = data_root |> parseABF
createDatasheet(all_files; filename=datafile, verbose=true)
dataset = openDatasheet(datafile, sheetName="all")
#updateDatasheet(datafile, all_files)
resB = ePhys.run_B_wave_analysis(dataset["All_Files"])
runAnalysis(datafile)
