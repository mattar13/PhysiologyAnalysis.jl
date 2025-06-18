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
main_channel = :red
img_fn3 = raw"G:\Data\Two Photon\2025-05-15-GRAB-DA_STR\b5_grabda-nircat-300uA_pulse014.tif"
stim_fn3 = raw"G:\Data\Patching\2025-05-15-GRAB-DA-STR\25515021.abf"

data_img = readImage(img_fn3)
stimulus = readABF(stim_fn3)

#%% We want to store the ROIs in a object
pixel_splits_roi!(data_img, 10)
make_circle_roi!(data_img, 1, 100, 100, 10)
getROIarr(data_img, 1)

#make a ROI that is a circle
#Pull out the framerate of the data
getSampleFreq(data_img)
