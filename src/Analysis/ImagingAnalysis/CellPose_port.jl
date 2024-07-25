using Revise
using Dates
using ElectroPhysiology 
using PhysiologyAnalysis
using Pkg; Pkg.activate("test")

using PyCall
# For some reason we need to do this weird dance. Just pyimport doesn't work and just python doesn't either
path_loc = splitpath(pathof(PhysiologyAnalysis))[1:end-1] |> joinpath
#ENV["CELLPOSE_LOCAL_MODELS_PATH"] = "C:/Users/mtarc/.julia/dev/PhysiologyAnalysis/src/Analysis/ImagingAnalysis/CellPoseModels"
py"""
import os
os.environ["CELLPOSE_LOCAL_MODELS_PATH"] = "C:/Users/mtarc/.julia/dev/PhysiologyAnalysis/src/Analysis/ImagingAnalysis/CellPoseModels"
import cellpose
from cellpose import models
"""
cellpose = pyimport("cellpose")
model = cellpose.models.Cellpose(model_type="cyto")

model