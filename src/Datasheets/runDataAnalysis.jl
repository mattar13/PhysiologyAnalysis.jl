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
datafile = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Members Folders\Brittney\JGP_data_analysis.xlsx"
all_files = data_root |> parseABF
#createDatasheet(all_files; filename=datafile, verbose=true)
#%% Make the changes talked about with paul
dataset = openDatasheet(datafile, sheetName="all")
#1) exclude all files with 0.5
dataset["All_Files"] = dataset["All_Files"] |> @filter(_.ND != 0.5) |> DataFrame 
#2) exclude all files with either a 5ms stimulus or a 16ms stimulus
dataset["All_Files"] = dataset["All_Files"] |> @filter(_.Stim_Time != 5 && _.Stim_Time != 16) |> DataFrame
#3) Change all green cone files to ND1 
for (idx, row) in enumerate(eachrow(dataset["All_Files"]))
     if row.Photoreceptor == "Cones" && row.ND == 0
          println(row)
          ND = 1
          WAVE = row.Wavelength
          PERCENT = row.Percent
          STIM_TIME = row.Stim_Time
          PHOTONS = ePhys.photon_lookup(WAVE,ND,PERCENT,ePhys.calibration_file) .* STIM_TIME
          println(PHOTONS)
          dataset["All_Files"][idx, :ND] = 1
          dataset["All_Files"][idx, :Photons] = PHOTONS
     end
end
dataset["All_Files"]
#%% 4) Include only those with certain values
dataset["All_Files"] |> @filter(_.Year == 2019) |> DataFrame

#%% Once you are ready save the XLSX file
XLSX.openxlsx(datafile, mode="w") do xf
     XLSX.rename!(xf[1], "All_Files") #Rename sheet 1
     try
          sheet = xf["All_Files"] #Try opening the All_Files
     catch
          println("Adding sheets")
          XLSX.addsheet!(xf, "All_Files")
     end
     XLSX.writetable!(xf["All_Files"],
          collect(DataFrames.eachcol(dataframe)),
          DataFrames.names(dataframe))
end


#updateDatasheet(datafile, all_files)
resB = ePhys.run_A_wave_analysis(dataset["All_Files"], verbose=true, 
     t_pre = 0.25, t_post = 0.5,
     a_cond="UNKNOWN", measure_minima = true
)
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