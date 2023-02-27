#%% Some common packages to use
using Revise, ePhys
using Pluto
ePhys.__init__()
run_experiment_analysis()
run_datasheet_analysis()
#Pluto.run()

#%% How can we make something that will remove lines in the 
#Open the dataframe
using Revise, ePhys
using DataFrames, Query, XLSX
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\Retinoschisis\\data_analysis.xlsx"
dataset = openDataset(datafile)
ePhys.__init__()
test = dataset["EXPERIMENTS"] |> @orderby(_.RSQ_fit) |> DataFrame
res = matchExperiment(dataset["TRACES"], test[1:10,:]) 

#%% Open Pauls files
ePhys.__init__()
using DataFrames, Query, XLSX
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)

#%% Try to save ABF
using Revise
using ePhys
import ePhys.saveABF
import ePhys.Experiment
using PyPlot
Revise.track(ePhys, "src/Readers/ABFReader/ABFReader.jl")
file_open = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Cones\2019_07_23_WT_P14_m1\Cones\Drugs\Green\nd0.5_1p_1ms\19723190.abf"
file_save = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Cones\2019_07_23_WT_P14_m1\Cones\Drugs\Green\nd0.5_1p_1ms\test.abf"
data = readABF(file_open, channels=-1) #This is necessary for saving
saveABF(data, file_save)