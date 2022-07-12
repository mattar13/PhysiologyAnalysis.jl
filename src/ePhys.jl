module ePhys

#println("Testing package works")
#=================== Here are the imports from other files ===================#
#using LsqFit #Used for fitting amplification and Intensity Response models
#using DSP
#using ContinuousWavelets
#using Wavelets
#using FFTW #Used for filtering
using Requires #This will help us load only the things we need
using Dates
using Base: String, println
import RCall as R #This allows us to use some R functionality
import PyCall as py #This allows us to use Python to call somethings (for )

export R, py

#=======================Import all experiment objects=======================#
include("Experiment/StimulusProtocol.jl")
include("Experiment/experiments.jl") #This file contains the Experiment structure. 

exp_readme = """
Experiment contains all of the necessary things to define an experiment. 

This must be included first before the files can be opened as it contains the base experiment object

"""

#======================Import all ABF extension imports======================#
include("Readers/ABFReader/ABFReader.jl") 
export readABF
export parseABF

#===============================ABF utilities===============================#
include("Utilities/ExperimentUtilities.jl")
include("Utilities/DataUtilities.jl")
export truncate_data, truncate_data!

#=Add filtering capability=#
using DSP
import Polynomials as PN #Import this (there are a few functions that get in the way)
include("Filtering/filtering.jl")
include("Filtering/filteringPipelines.jl")
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

#===============================Import all Datasheet tools==============================#
using DataFrames, Query, XLSX
include("Datasheets/DatasheetFunctions.jl")

#====================Import all the tools needed to analyze the data====================#
#First import models necessary for the analysis

using Statistics, StatsBase #These functions use R functions as well as StatsBase
include("Analysis/Stats.jl")
export RSQ

using Distributions
include("Analysis/Models.jl")
include("Analysis/ERGAnalysis.jl")
export calculate_basic_stats
export saturated_response, dim_response
export minima_to_peak, time_to_peak
export percent_recovery_interval #This finds the dominant time constant
export recovery_time_constant #This finds the recovery time constant
export integral #This finds the integration time
export get_response
export amplification
export curve_fit #curve fitting from LsqFit
export IR_curve

include("Analysis/TimescaleAnalysis.jl")
export calculate_threshold
export get_timestamps
export max_interval_algorithim
export timescale_analysis

#========================================Plotting utilities========================================#
include("Plotting/PhysPlotting.jl")
export plot, plot!
export plt

end