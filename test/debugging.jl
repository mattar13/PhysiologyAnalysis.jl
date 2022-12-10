using Revise
using ePhys
using PyPlot
using DataFrames, Query, XLSX
ePhys.__init__()

#%% Why is the data from wildtype looking so weird
#We only want to filter one channel
file_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_04_21_a13MelCreAdult\Mouse2_Adult_WT\BaCl_LAP4\Rods"
files = file_root |> parseABF
@time data = data_filter(readABF(files), channels = -1);

fig, ax = plt.subplots(2)
plot_experiment(ax[1], data, channels = 1, c = :black)
plot_experiment(ax[2], data, channels = 2, c = :black)

#Split the channels
channels = ePhys.eachchannel(data) |> collect

#%% calculate rmaxes
rmaxes = saturated_response(dataBASE)
p_rec = percent_recovery_interval(dataBASE, rmaxes)
maximum(p_rec)
#%% Open Pauls files
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)