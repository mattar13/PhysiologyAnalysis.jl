using Revise
using ePhys
using PyPlot
using DataFrames, Query, XLSX
ePhys.__init__()

# Why is the data from wildtype looking so weird
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2021_09_22_C59S-13\Mouse3_P13_C59S\BaCl_LAP4\Rods\nd4_1p_0000.abf"
root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_07_11_Adult\Mouse1_Adult_WT"
files_WT30a = "$root\\BaCl_LAP4\\Rods" |> parseABF
@time data = data_filter(readABF(files_WT30a, channels = ["Vm_prime"]));
data.stim_protocol 


#%% 
data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Projects\Retinoschisis" #The data root
data_file = joinpath(data_root, "data_analysis.xlsx") #This is the main file we can use. The root may change
dataset = openDatasheet(data_file, sheetName="all")
RS_trace_B = dataset["trace_B"]
P30_WT = (2022, 7, 11, 1, "Vm_prime4") #P30 WT (2021_9_28_n2_Vm_prime_Rods)
q_WT30b = matchExperiment(RS_trace_B, P30_WT) |> DataFrame
t_post_val = 2.0
@time data_WT30a, data_WT30ab = data_filter.(readABF(q_WT30b), t_post=t_post_val)
data_WT30B = data_WT30ab - data_WT30a

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

