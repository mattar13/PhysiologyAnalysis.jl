using Revise
using ePhys
using PyPlot
using DataFrames, Query, XLSX
ePhys.__init__()
import ePhys.CWTprocess
using ContinuousWavelets
using Statistics
#%% enter in the file path of the file you would like to analyze here
path = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_12_20_Matt_WTAdult\Mouse2_WT_Adult\BaCl_LAP4\Rods\nd3_1p_0002.abf"
channels = ["Vm_prime"]
data = readABF(path, channels=channels) |> average_sweeps
baseline_adjust!(data)
ePhys.truncate_data!(data, t_pre=0.1, t_post=4.0)
data = data * -1
# Look more into how the wavelets could be used
wave = ePhys.Morlet(0.50π)
period_window = (2^0, 2^4)
power_window = (0.0, 1.0)
dataRE, CWT_OG = cwt_filter(data,
     wave=wave,
);

dataARTIFACT, CWT_ART = cwt_filter(data,
     wave=wave,
     period_window=period_window, power_window=power_window
);


CWTo_PROC = CWTprocess(CWT_OG[1, :, :, 1])
mu = mean(CWTo_PROC[.!isinf.(CWTo_PROC)])
sig = std(CWTo_PROC[.!isinf.(CWTo_PROC)])
levelsO = LinRange(mu-2*sig, mu+2*sig, 40)
freqsO = log.(2, 1:size(CWT_OG, 3))

CWTa_PROC = CWTprocess(CWT_ART[1, :, :, 1])
mu = mean(CWTa_PROC[.!isinf.(CWTa_PROC)])
sig = std(CWTa_PROC[.!isinf.(CWTa_PROC)])
levelsA = LinRange(mu - 2*sig, mu + 2*sig, 40)
freqsA = log.(2, 1:size(CWT_ART, 3))

#%% Plot the results
plt.close("all")
figTEST, axTEST = plt.subplots(2, 2, figsize=(15, 7.4))
imF1 = axTEST[1, 1].contourf(data.t, freqsO, CWTo_PROC,
     cmap="seismic", levels=levelsA, extend = "both", 
)
cbar_ax1 = figTEST.add_axes([0.925, 0.5, 0.03, 0.25])
figTEST.colorbar(imF1, cax=cbar_ax1, orientation="vertical")
plot_experiment(axTEST[2, 1], data * -1, channels=1, c=:black)
plot_experiment(axTEST[2, 1], dataRE, channels=1, c=:red)

imF2 = axTEST[1, 2].contourf(data.t, freqsA, CWTa_PROC,
     cmap="seismic", levels=levelsA, extend = "both", 
)
cbar_ax2 = figTEST.add_axes([0.025, 0.5, 0.03, 0.25])
figTEST.colorbar(imF2, cax=cbar_ax2, orientation="vertical")
plot_experiment(axTEST[2, 2], data * -1, channels=1, c=:black)
plot_experiment(axTEST[2, 2], dataARTIFACT * -1, channels=1, c=:red)

#%% Open Pauls files
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)