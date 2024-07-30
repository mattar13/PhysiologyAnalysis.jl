"""
In order to run this, we need to have python included in the path
"""
function build_cellpose(;python_env = "")
     ENV["PYTHON"] = python_env
     Pkg.build("PyCall")
     #This is necessary to run cellpose from new
     Conda.pip_interop(true)
     Conda.pip("install", "numpy")
     Conda.pip("install --user", "opencv-python")
     Conda.pip("install", "cellpose")
end

function cellpose_model(;model_type="cyto", relative_path_loc = "Analysis\\ImagingAnalysis\\CellPoseModels\\")
     #╔═╡Set up the python environment to play nice with julia
     path_loc = joinpath(splitpath(pathof(PhysiologyAnalysis))[1:end-1]..., relative_path_loc) 
     try
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
     catch error
          #If the error is with python, I want to build_cellpose from PyCall
          throw(error)
     end
end
