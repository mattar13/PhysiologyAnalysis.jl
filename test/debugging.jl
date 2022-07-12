using Revise
using PhysAnalysis
plt.pygui(true)

#%%Lets try some file opening and analysis
folder = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Organoids\2022_06_17_Organoid"
files = parseABF(folder)
#%%
data = readABF(files[1])
average_sweeps!(data)
baseline_adjust!(data; mode = :mean) #This centers the data
baseline_adjust!(data; polyN = 1, region = :whole) #This removes drift
truncate_data!(data, t_post = 1.0)
lowpass_filter!(data, freq = 300.0)
#highpass_filter!(data, freq = 0.05)
#EI_filter!(data)

size(data)
plot(fig1, data, ylims = (-0.005, 0.005))

#%% Open Pauls files
using MAT
file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\MAT files\2022-Feb-26_RBC_SPR.mat"
data = matopen(file)
vars = matread(file)