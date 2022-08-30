using Revise
using ePhys
using DataFrames
#%% Using the dataframe utilites to make Datasheets

file_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Noise and ArtifactTests\2022_08_29_VolumeToRMSNoise"
noise_old_chamber = "$(file_root)/I0 chamber_0000.abf"
noise_2_4 = "$(file_root)/I0 2.4ml 3.15ID.abf"
data = readABF(noise_2_4, channels = ["Vm_prime"])


#To run the analysis on a single dataset this is the max code that is needed
data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis" #The data root
datafile = "$(data_root)/data_analysis.xlsx"
all_files = data_root |> parseABF
dataset = openDatasheet(datafile, sheetName = "all")

df, resA, resB, resG = runAnalysis(datafile, all_files)

#%% Test plotting
using PyPlot
import PyPlot.plt
pygui(true) #Activate the GUI so the plots are not inline
import ePhys: truncate_data!, average_sweeps!, baseline_adjust!
import ePhys.fft_spectrum
import ePhys: filter_data, dwt_filter


#%% Open Pauls files
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)

