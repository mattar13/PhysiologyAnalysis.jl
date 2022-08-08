using Revise
using ePhys
using DataFrames
#%% Test plotting
using PyPlot
import PyPlot.plt
import ePhys: Lowpass, Highpass, Butterworth, digitalfilter, filt
pygui(true) #Activate the GUI so the plots are not inline


#%% Lets try to make a Pluto interactive notebook that will help analyze the data




#%% Using the dataframe utilites to make Datasheets
ePhys.__init__()
data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis" #The data root
all_files = data_root |> parseABF

datasheet = ePhys.createDatasheet(all_files)
ePhys.run_A_wave_analysis(datasheet)
ePhys.run_B_wave_analysis(datasheet)
ePhys.run_G_wave_analysis(datasheet)

#%% Open Pauls files
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)