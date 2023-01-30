using Revise
using ePhys
using Pluto

#Pluto.run(notebook = "src/Interface/test_interface.jl")

#%% This section will determine filter settings
Pluto.run(notebook="src/Interface/experiment_analysis.jl")
Pluto.run(notebook="src/Interface/file_analysis.jl") #A-wave analysis for a file of traces
Pluto.run(notebook = "src/Interface/filter_determination.jl") #A-wave analysis for a file of traces
Pluto.run(notebook = "src/Interface/subtraction_analysis.jl")