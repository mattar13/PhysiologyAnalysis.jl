#%% Section 1. Revision of interface
using Revise, PhysiologyAnalysis
using Pkg; Pkg.activate("./test/")

#%%
using Pluto
run_experiment_analysis()
run_datasheet_analysis()
#Pluto.run()

#%% Section 2, Revisions of matching experiments
#Open the dataframe
using Revise, PhysiologyAnalysis
using DataFrames, Query, XLSX
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\Retinoschisis\\data_analysis.xlsx"
dataset = openDataset(datafile)
PhysiologyAnalysis.__init__()
test = dataset["EXPERIMENTS"] |> @orderby(_.RSQ_fit) |> DataFrame
res = matchExperiment(dataset["TRACES"], test[1:10,:]) 

#%% Section 3, Opening Matlab IRIS files
PhysiologyAnalysis.__init__()
using DataFrames, Query, XLSX
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)

#%% Section 4, Saving ABF files
using Revise
using PhysiologyAnalysis
import PhysiologyAnalysis.saveABF
import PhysiologyAnalysis.Experiment
using PyPlot
Revise.track(PhysiologyAnalysis, "src/Readers/ABFReader/ABFReader.jl")
file_open = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Cones\2019_07_23_WT_P14_m1\Cones\Drugs\Green\nd0.5_1p_1ms\19723190.abf"
file_save = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Cones\2019_07_23_WT_P14_m1\Cones\Drugs\Green\nd0.5_1p_1ms\test.abf"
data = readABF(file_open, channels=-1) #This is necessary for saving
saveABF(data, file_save)