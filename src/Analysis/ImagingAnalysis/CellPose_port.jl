using PyCall, Conda

function build_cellpose(;python_env = "", model_dir = nothing)
     if isnothing(model_dir)
          error("Please specify a directory for the models to reside in")
     else
          ENV["PYTHON"] = python_env
          Pkg.build("PyCall")
          #This is necessary to run cellpose from new
          Conda.pip_interop(true)
          Conda.pip("install", "opencv-python")
          Conda.pip("install")
     end
end

function cellpose_model(;model_type="cyto", relative_path_loc = "Analysis\\ImagingAnalysis\\CellPoseModels\\")
     #╔═╡Set up the python environment to play nice with julia
     path_loc = joinpath(splitpath(pathof(PhysiologyAnalysis))[1:end-1]..., relative_path_loc) 
     py"""
     import os
     os.environ["CELLPOSE_LOCAL_MODELS_PATH"] = $path_loc
     import cellpose
     from cellpose import models
     """
     #╔═╡Import and create the models
     cellpose = pyimport("cellpose")
     model = cellpose.models.Cellpose(model_type=model_type)
     return model
end

function cellpose_eval(img)
     mask, flow, style, diam = model.eval(grn_zproj)
     return mask
end

