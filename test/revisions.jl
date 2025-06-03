using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
using Pkg; Pkg.activate("test")
using GLMakie, PhysiologyPlotting

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
println("Loading the quinpirole baseline data...")
img_da_10mM_fn = raw"H:\Data\Two Photon\2025-05-19-nirCAT-STR-DPUFF\nirCAT_s2_10mM_mq010.tif"
stim_da_10mM_fn = raw"H:\Data\Patching\2025-05-19-nirCAT-str\25519016.abf"

experiment = readImage(img_da_10mM_fn)
stim_trace = readABF(stim_da_10mM_fn, flatten_episodic = true,stimulus_name = "IN 2", stimulus_threshold = 0.5)
stimulusProtocol = getStimulusProtocol(stim_trace)
addStimulus!(experiment, stimulusProtocol)
stim_ends = getStimulusEndTime(experiment)
dt = experiment.dt

