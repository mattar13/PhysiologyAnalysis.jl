using Revise, ePhys
using DataFrames, Query, XLSX
using Dates
import ePhys: openDatasheet, run_B_wave_analysis
yyyy = year(Dates.now())
mm = month(Dates.now())
dd = day(Dates.now())

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
data_root = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\JGP_Files\Traces"
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\GNAT\\$(yyyy)_$(mm)_$(dd)_cone_data_analysis.xlsx"
all_files = data_root |> parseABF
#createDatasheet(all_files; filename=datafile, verbose=true)
dataset = openDatasheet(datafile, sheetName="all")
#updateDatasheet(datafile, all_files)
resB = ePhys.run_A_wave_analysis(dataset["All_Files"], verbose=true, 
     t_pre = 0.25, t_post = 0.5,
     a_cond="UNKNOWN", measure_minima = true
)
datafile = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Members Folders\Brittney\JGP_data_analysis.xlsx"
ePhys.add_analysis_sheets(resB, datafile; append = "B")
runAnalysis(datafile)

#%% Run some other files
data_files1 = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Cones" |> parseABF
data_files2 = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Rods" |> parseABF
data_files3 = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\JGP_Files\Traces" |> parseABF
data_files4 = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Restructured" |> parseABF
all_files = [data_files1..., data_files2..., data_files4...]
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\GNAT\\cone_data_analysis.xlsx"
updateDatasheet(datafile, all_files)
runAnalysis(datafile)
#%%

test_file