#%% Using the Julia interface
using Revise, ePhys 
using Pluto
run_experiment_analysis()

#%% Something wrong with the experiment analysis and models
using Revise, ePhys
using DataFrames, Query, XLSX
import ePhys: run_A_wave_analysis, run_B_wave_analysis, run_G_wave_analysis
#enter in the file path of the file you would like to analyze here
exp_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_07_11_Adult\Mouse3_Adult_RS1KO"
experiment_paths = exp_root |> parseABF
savefile_name = nothing
all_files = createDatasheet(experiment_paths; filename = savefile_name)
all_files
res_A = run_A_wave_analysis(all_files)
res_B = run_B_wave_analysis(all_files)
res_G = run_G_wave_analysis(all_files)

#%% Here we want to test out the convienance functions
using DataFrames, Query, XLSX
import ePhys: GenerateFitFrame, STFfit

data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Projects\GNAT"
data_file = joinpath(data_root, "data_analysis.xlsx") #This is the main file we can use. The root may change
GN_trace_B = openDatasheet(data_file; sheetName="trace_B")

p14_TRACE = GN_trace_B |> @filter(_.Age == 14  && _.Wavelength == 520) |> DataFrame#Extract all data for P14
lb_STF = (100.0, 1.0, 0.1)
p0_STF = (500.0, 200.0, 2.0)
ub_STF = (3000.0, 3000.0, 10.0)
p14_STF = GenerateFitFrame(p14_TRACE, :Response, :Maxima; lb = lb_STF, p0 = p0_STF, ub = ub_STF) #|> @filter(_.RSQ > 0.50) |> @orderby(_.RSQ)  |> DataFrame
#(rmax, k, n)
lb_IR = (1.0,    10^-1, 0.1 )
p0_IR = (500.0,  200.0, 2.0 )
ub_IR = (2400.0, 10^6,  10.0)
p14_IR = GenerateFitFrame(p14_TRACE, :Photons, :Response; lb = lb_IR, p0 = p0_IR, ub = ub_IR) #|> @filter(_.RSQ > 0.50) |> @orderby(_.RSQ)  |> DataFrame

p14_TRACE.Photons

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

#%% Open Pauls files
ePhys.__init__()
using DataFrames, Query, XLSX
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)