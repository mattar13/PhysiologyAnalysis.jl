using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON
# using Pkg; Pkg.activate("test")
# using GLMakie, PhysiologyPlotting

#%% ╔═╡Found an issue. With really crazy spikes, the baseline correction is not working.
println("Loading the quinpirole baseline data...")
img_fn = raw"F:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b6_grabda-nircat-300uA_pulse016.tif"
stim_fn = raw"F:\Data\Patching\2025-05-15-GRAB-DA-STR\25515025.abf"

data_quin_base = load_electric_data(img_fn, stim_fn, pre_event_time = 20.0, post_event_time = 60.0, main_channel = :grn)