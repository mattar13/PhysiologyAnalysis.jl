module PhysAnalysis

println("Testing package works")
#=================== Here are the imports from other files ===================#
#using Statistics
#using Polynomials
#using Distributions
#using StatsBase #Used for polynomial fitting
#using LsqFit #Used for fitting amplification and Intensity Response models
#using DSP
#using ContinuousWavelets
#using Wavelets
#using FFTW #Used for filtering

#======================Import all ABF extension imports======================#
using ABFReader

#=======================Import all experiment objects=======================#
include("Experiment/StimulusProtocol.jl")
include("Experiment/Experiments.jl") #This file contains the Experiment structure. 


#===============================ABF utilities===============================#
include("OpeningFiles/OpeningABF.jl")
export readABF


include("Utilities/ExperimentUtilities.jl")
include("Utilities/DataUtilities.jl")
export truncate_data, truncate_data!

#=Add filtering capability=#
using Polynomials, DSP
include("Analysis/Filtering.jl")
#export filter_data #Don't export this one explicitly
export baseline_adjust, baseline_adjust!
export lowpass_filter, lowpass_filter!
export highpass_filter, highpass_filter!
export notch_filter, notch_filter!
export EI_filter, EI_filter!
export cwt_filter, cwt_filter!
export dwt_filter
export average_sweeps, average_sweeps!
export normalize, normalize!

#====================Import all the tools needed to analyze the data====================#
using Distributions, StatsBase
include("Analysis/ERGAnalysis.jl")
export RSQ
export calculate_basic_stats
export saturated_response, dim_response, minima_to_peak, time_to_peak
export get_response
export pepperburg_analysis
export integral
export recovery_tau
export amplification
export curve_fit #curve fitting from LsqFit
export IR_curve

include("Analysis/WholeCellAnalysis.jl")
export calculate_threshold
export get_timestamps
export max_interval_algorithim
export timescale_analysis

end