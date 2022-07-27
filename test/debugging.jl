using Revise
using ePhys
using JLD2
import ePhys.plt: xlim, ylim
import ePhys.plot_experiment

#Test the timeseries analysis used for modelling
dt = 0.5
file_loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Patching"
target_file = "$(file_loc)/2019_11_03_Patch/Animal_2/Cell_3/19n03042.abf"
#Eventually we will open multiple analysis files
data = readABF(target_file, channels=["Vm_prime4"], stimulus_name=nothing, time_unit=:ms)
data.t
down_data = downsample(data, 1/dt)
down_data.t
down_data.data_array
downsample!(data, 1/0.5)
data.t
#add a filter setting for downsample

region = (-100, 2500)
data_root = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/figure_data"
isolated_path = "$(data_root)/isolated_model"
isolated_data = load("$(isolated_path)/data.jld2")
isolated_timestamps = load("$(isolated_path)/timestamps.jld2")

timeseries_analysis(isolated_data["Time"], isolated_data["DataArray"])

#%% Make the figure here reqlly quick
#%% Lets pull up some good wildtype data
abg_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_07_11_Adult\Mouse2_Adult_Abino\NoDrugs"
ab_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_07_11_Adult\Mouse2_Adult_Abino\BaCl"
a_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_07_11_Adult\Mouse2_Adult_Abino\BaCl_LAP4"

abg_path = parseABF(abg_root)
ab_path = parseABF(ab_root)[2:end]
a_path = parseABF(a_root)
import ePhys.data_filter

abg_data = readABF(abg_path) |> data_filter
ab_data = readABF(ab_path) |> data_filter
a_data = readABF(a_path) |> data_filter
b_data = ab_data - a_data
g_data = abg_data - ab_data
#%%
width_inches = 10.0
height_inches = 5.0
fig2 = plt.figure("Methods", figsize=(width_inches, height_inches))
#fig2.text(0.01, 0.90, "A", va="center", fontsize=20.0)
#fig2.text(0.01, 0.40, "B", va="center", fontsize=20.0)
gs1 = fig2.add_gridspec(2, 1,
     right=0.95, left=0.05,
     top=0.93, bottom=0.0,
     wspace=1.0, hspace=0.40
) #Make a 3x3 figure
#fig2.tight_layout(pad=5.0)

ylims_i = (-1000, 2000) #use these as the fig2 ylims_i
xlims_i = (-0.3, 1.5) #use this for the fig2 xlims_i

fig2_A = gs1[1, 1].subgridspec(ncols=3, nrows=1) #row 1 is 3 spots
axA1 = fig2.add_subplot(fig2_A[1, 1])
xlim(xlims_i);
ylim(ylims_i);
#title("Ringers Solution", fontsize=20.0)
plot_experiment(axA1, abg_data, c=:black, lw=2.0)
axA1.vlines(0.0, ymin=ylims_i[1], ymax=ylims_i[2], linestyle="dashed", colors=:black, lw=2.0)
#add_scalebar(fig2_A1, (-0.2, -700.0), (0.4, 250.0), fontsize=8.0)

axA2 = fig2.add_subplot(fig2_A[1, 2])
xlim(xlims_i);
ylim(ylims_i);
#title("Blocking Glial Cells", fontsize=20.0)
plot_experiment(axA2, ab_data, c=:black, lw=3.0)
axA2.vlines(0.0, ymin=ylims_i[1], ymax=ylims_i[2], linestyle="dashed", colors=:black, lw=2.0)
#fig2_A2.annotate("+ BaCl\$_2\$", [0.5, 250.0], va="center", ha="center", fontsize=10.0)

axA3 = fig2.add_subplot(fig2_A[1, 3])
xlim(xlims_i);
ylim(ylims_i);
#title("Blocking Bipolar Cells", fontsize=20.0)
plot_experiment(axA3, a_data, c=:black, lw=3.0)
fig2_A3.vlines(0.0, ymin=ylims_i[1], ymax=ylims_i[2], linestyle="dashed", colors=:black, lw=2.0)

#fig2_A3.annotate("+ LAP-4, BaCl\$_2\$, Asp", [0.8, 250.0], va="center", ha="center", fontsize=10.0)

# Row B (two plots are offset)
fig2_B = gs1[2].subgridspec(ncols=4, nrows=1, width_ratios=[0.1, 0.4, 0.4, 0.1]) #row 2 is 2

fig2_B1 = fig2.add_subplot(fig2_B[1, 2]) #change the spacing on this one
xlim(xlims_i);
ylim(ylims_i);
plot_experiment(fig2_B1, g_data, color=:black, lw=3.0)
fig2_B1.vlines(0.0, ymin=-1000.0, ymax=800.0, linestyle="dashed", colors=:black, lw=2.0)
#fig2_B1.annotate("Isolated Glial Component", [0.8, 100.0], va="center", ha="center", fontsize=10.0)

fig2_B2 = fig2.add_subplot(fig2_B[1, 3])
xlim(xlims_i);
ylim(ylims_i);
plot_experiment(fig2_B2, b_data, color=:black, lw=3.0)
fig2_B2.vlines(0.0, ymin=-1000.0, ymax=800.0, linestyle="dashed", colors=:black, lw=2.0)

#%%Some coding to run analysis
#lets test conditional loading of Pyplot files
@require Plot begin
     println("GR plot is loaded")
end

#%%Lets try some file opening and analysis
folder = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Organoids\2022_06_17_Organoid"
files = parseABF(folder)
#%%
data = readABF(files[1])
average_sweeps!(data)
baseline_adjust!(data; mode=:mean) #This centers the data
baseline_adjust!(data; polyN=1, region=:whole) #This removes drift
truncate_data!(data, t_post=1.0)
lowpass_filter!(data, freq=300.0)
#highpass_filter!(data, freq = 0.05)
#EI_filter!(data)

size(data)
plot(fig1, data, ylims=(-0.005, 0.005))

#%% Open Pauls files
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)