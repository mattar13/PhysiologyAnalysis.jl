module PhysiologyAnalysis

using Requires
# The top level is the ElectroPhysiology package. These are not imported into the workspace
using Dates

using ElectroPhysiology
import ElectroPhysiology: Experiment, readABF, parseABF
import ElectroPhysiology: now, year, month, day, hour, minute, second

#= Packages used for fitting data ====================================#
using LsqFit #Used for fitting amplification, Intensity Response, and Resistance Capacitance models

#= Packages used for Analyzing data ==================================#
import Polynomials as PN #used for fitting and stats
#@time using DataFrames, Query, XLSX #Load these extra utilites immediately
import XLSX: readtable, readxlsx #Import XLSX commands

using StatsBase #Used for mean, std, and SEM functions.
using HypothesisTests

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

#3) Import all Datasheet tools ===========================================================#
function __init__()
     @require OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed" begin
          @require ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78" begin
               println("Loading differential equations and nose fitting")
               using ModelingToolkit
               using OrdinaryDiffEq
               include("Fitting/NoseModel.jl")
               export findNosePeak
          end
     end

     @require DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0" begin
          using DataFrames
          using Query, XLSX
          
          println(Query)
          
          println(XLSX)
          
          export readtable, readxlsx, XLSX
          include("Datasheets/RegexFunctions.jl")

          include("Datasheets/FilePathExtraction.jl")

          include("Datasheets/DatasheetFunctions.jl")
          export matchDataset, excludeDataset, concatDatasets

          include("Datasheets/DatasheetCreation.jl")
          export openDataset, createDataset
          export saveDataset, backupDataset

          include("Datasheets/DatasheetAnalysis.jl")
          export runAnalysis
          export runTraceAnalysis
          export runExperimentAnalysis
          export runConditionsAnalysis
          export runStatsAnalysis
          export matchExperiment, excludeExperiment
          export flagExperiment, flagExperiment!, unflagALL!
          export parseColumn!
          export GenerateFitFrame

          include("Datasheets/extra_utilities.jl")
     end
end

end