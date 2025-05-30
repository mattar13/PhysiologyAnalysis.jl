using Dates
using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics, NaNMath
import ElectroPhysiology: Experiment, TWO_PHOTON
using Pkg; Pkg.activate("test")
using GLMakie

#%% ╔═╡This task is for extraction of points, centroids, and ROIs using cellpose
println("Loading the quinpirole baseline data...")
img_fn = raw"F:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b6_grabda-nircat-300uA_pulse016.tif"
stim_fn = raw"F:\Data\Patching\2025-05-15-GRAB-DA-STR\25515025.abf"

data_quin_base = load_electric_data(img_fn, stim_fn, red_lam = 1e3, grn_lam = 1e3)
all_rois = data_quin_base["mean_sig_trace"]
