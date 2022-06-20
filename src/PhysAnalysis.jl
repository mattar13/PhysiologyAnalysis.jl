module PhysAnalysis

println("Testing package works")

#=======================Import all experiment objects=======================#
include("Analysis\\Experiments.jl") #This file contains the Experiment structure. 

include("Analysis\\filtering.jl")
#export filter_data #Don't export this one explicitly
export baseline_cancel, baseline_cancel!
export lowpass_filter, lowpass_filter!
export highpass_filter, highpass_filter!
export notch_filter, notch_filter!
export EI_filter, EI_filter!
export cwt_filter, cwt_filter!
export dwt_filter
export average_sweeps, average_sweeps!
export normalize, normalize!

#======================Import all ABF extension imports======================#
using ABFReader
export readABF #Eventually this will be folded into a file that automatically determines how to open the extension


end