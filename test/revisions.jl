using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
println("Loading the quinpirole baseline data...")
img_da_10mM_fn = raw"H:\Data\Two Photon\2025-05-19-nirCAT-STR-DPUFF\nirCAT_s2_10mM_mq010.tif"
stim_da_10mM_fn = raw"H:\Data\Patching\2025-05-19-nirCAT-str\25519016.abf"
data_da_10mM = load_puffing_data(img_da_10mM_fn, stim_da_10mM_fn, split_channel = false, main_channel = :grn)



mean_sig_traces_da_10mM = mean(data_da_10mM["sig_traces"], dims = 2)[:,1,:,:]
mean_sig_trace_da_10mM = data_da_10mM["mean_sig_trace"]
time_da_10mM = data_da_10mM["time"][1:size(mean_sig_trace_da_10mM, 1)]


#Marla wants me to do a different baseline analysis. Lets work on some new ways to do this. Fitting an exponential to the baseline.


