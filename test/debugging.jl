using Revise
using PhysAnalysis

#file to open
test_file = "test\\to_analyze.abf"
data = readABFInfo(test_file)


target_paths = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Organoids\2022_06_17_Organoid"
paths = target_paths |> parse_abf;


#%% Fix the patch clamp measures
target_folder = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Patching\Jordans_Patch_Data\UsuableData"
paths = target_folder |> parse_abf;

path = paths[1]
data = readABF(path, channels=["Im_scaled"], stimulus_name=nothing, flatten_episodic=true)
threshold = calculate_threshold(data, Z=4)[:, 1, 1]
spike_timestampsi = ABFReader.get_timestamps(data)

#%%
burst_timestampsi = ABFReader.max_interval_algorithim(spike_timestampsi, verbose=true)

#%%
spike_durs, isi = ABFReader.extract_interval(spike_timestampsi)
burst_durs, ibi = ABFReader.extract_interval(burst_timestampsi[1])

burst_timestampsi

plot(data.t[t_idxs[1]:100:t_idxs[2]], data.data_array[:, t_idxs[1]:100:t_idxs[2], 1]')
hline!([threshold[1, 1, 1]])
vline!(spike_t_start, c=:green)
vline!(spike_t_end, c=:red)


spike_durs
spike_timestampsi[1]
#%% Open Pauls files
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)