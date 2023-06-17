module PhysiologyAnalysis

# The top level is the ElectroPhysiology package. These are not imported into the workspace
using Dates

using ElectroPhysiology
import ElectroPhysiology: Experiment, readABF, parseABF
import ElectroPhysiology: now, year, month, day, hour, minute, second

#= Packages used for fitting data ====================================#
using LsqFit #Used for fitting amplification, Intensity Response, and Resistance Capacitance models

#= Packages used for Analyzing data ==================================#
import Polynomials as PN #used for fitting and stats
using DataFrames, Query, XLSX #Load these extra utilites immediately
import XLSX: readtable, readxlsx #Import XLSX commands

using StatsBase #Used for mean, std, and SEM functions.
using HypothesisTests
using ModelingToolkit, OrdinaryDiffEq
#= Packages not yet uses
using Distributions
using Statistics, StatsBase #These functions use R functions as well as StatsBase

#export some basic functions from ============================================================#
export readABF, parseABF
export plt

#This package does 3 things: 

#1) Fitting ============================================================================#
include("Fitting/Models.jl")
export HILL_MODEL, HILLfit, STFfit
export AMP, AMPfit
export curve_fit #curve fitting from LsqFit

include("Fitting/NoseModel.jl")
export findNosePeak

#2) Data anlysis ========================================================================#
include("Analysis/ERGAnalysis.jl")
#export calculate_basic_stats
export saturated_response, dim_response
export minima_to_peak, time_to_peak
export percent_recovery_interval #This finds the dominant time constant
export recovery_time_constant #This finds the recovery time constant
export integral #This finds the integration time
export get_response
export calculate_threshold

#using JLD2 #Maybe this should be added to Requires.jl
include("Analysis/TimescaleAnalysis.jl")
export get_timestamps, extract_interval
export max_interval_algorithim, timeseries_analysis

include("Analysis/Stats.jl")
export dataset_statistics



#3) Import all Datasheet tools ===========================================================#
export readtable, readxlsx, XLSX
include("Datasheets/RegexFunctions.jl")

include("Datasheets/FilePathExtraction.jl")

include("Datasheets/DatasheetFunctions.jl")

include("Datasheets/DatasheetCreation.jl")
export openDataset, createDataset
export saveDataset, backupDataset

include("Datasheets/DatasheetAnalysis.jl")
export runAnalysis
export runTraceAnalysis
export runExperimentAnalysis
export runConditionsAnalysis
export matchExperiment, excludeExperiment
export flagExperiment, flagExperiment!, unflagALL!
export parseColumn!
export GenerateFitFrame

include("Datasheets/extra_utilities.jl")
#5) Plotting utilities will be loaded in automatically ==============================================#

#= This may be included in requires in PhysiologyPlotting
include("Plotting/DatasheetPlotting.jl")
export plot_ir_scatter, plot_ir_fit, plot_IR
export plot_data_summary
export plot_dataset_fits, plot_dataset_vals
=# 

end