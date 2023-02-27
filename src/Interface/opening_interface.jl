#Pluto.run(notebook = "src/Interface/test_interface.jl")
run_experiment_analysis() = Pluto.run(notebook="src/Interface/experiment_analysis.jl")
run_datasheet_analysis() = Pluto.run(notebook="src/Interface/datasheet_analysis.jl")
run_trace_analysis() = Pluto.run(notebook="src/Interface/file_analysis.jl") #A-wave analysis for a file of traces
run_filter_determination() = Pluto.run(notebook = "src/Interface/filter_determination.jl") #A-wave analysis for a file of traces
run_subtraction_analysis() = Pluto.run(notebook = "src/Interface/subtraction_analysis.jl")