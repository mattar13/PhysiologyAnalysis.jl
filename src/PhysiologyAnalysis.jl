module PhysiologyAnalysis

# The top level is the ElectroPhysiology package
using ElectroPhysiology
export readABF
export Experiment

#=================== Here are the imports from other files ===================#
using Requires #This will help us load only the things we need
using Crayons #Really cool package for coloring text for debugging

#=Add filtering capability=#
using DSP #Used for lowpass, highpass, EI, and notch filtering
using LsqFit #Used for fitting amplification, Intensity Response, and Resistance Capacitance models
import Polynomials as PN #Import this (there are a few functions that get in the way)
include("Filtering/filtering.jl")
#include("Filtering/filteringPipelines.jl") #Not ready to uncomment this one yet
#export filter_data #Don't export this one explicitly
export baseline_adjust, baseline_adjust!
export lowpass_filter, lowpass_filter!
export highpass_filter, highpass_filter!
export notch_filter, notch_filter!
export EI_filter, EI_filter!
export cwt_filter, cwt_filter!
export dwt_filter

#export rolling_mean
#export normalize, normalize!

include("Filtering/filteringPipelines.jl")
export data_filter!, data_filter

include("Fitting/fitting.jl")
export MeanSquaredError

#====================Import all the tools needed to analyze the data====================#
#First import models necessary for the analysis
using Statistics, StatsBase #These functions use R functions as well as StatsBase
include("Analysis/Stats.jl")
export RSQ

using Distributions
include("Analysis/Models.jl")
export HILL_MODEL
include("Analysis/ERGAnalysis.jl")
#export calculate_basic_stats
export saturated_response, dim_response
export minima_to_peak, time_to_peak
export percent_recovery_interval #This finds the dominant time constant
export recovery_time_constant #This finds the recovery time constant
export integral #This finds the integration time
export get_response
export amplification
export curve_fit #curve fitting from LsqFit
export IR_curve
export calculate_threshold

#using JLD2 #Maybe this should be added to Requires.jl
#include("Analysis/TimescaleAnalysis.jl")
#export get_timestamps, extract_interval
#export max_interval_algorithim, timeseries_analysis

#Dataframe utilities are baked in automatically
using DataFrames, Query, XLSX #Load these extra utilites immediately
import XLSX: readtable, readxlsx #Import XLSX commands
export readtable, readxlsx, XLSX
include("Datasheets/RegexFunctions.jl")
include("Datasheets/FilePathExtraction.jl")
include("Datasheets/DatasheetFunctions.jl")
include("Datasheets/DatasheetCreation.jl")
include("Datasheets/DatasheetAnalysis.jl")
export openDataset, createDataset, updateDataset
export runAnalysis
export runTraceAnalysis
export matchExperiment
export parseColumn!
export GenerateFitFrame
export saveDataset, backupDataset

#Plotting utilities will be loaded in automatically
import PyPlot
import PyPlot.plt #All the base utilities for plotting
import PyPlot.matplotlib
import PyCall as py #This allows us to use Python to call somethings 
import PyCall: @pyimport, PyObject

include("Plotting/DefaultSettings.jl") #This requires PyPlot
include("Plotting/PlottingUtilities.jl")
include("Plotting/PhysPyPlot.jl")
export plot_experiment, plot_experiment_fit


#using DataFrames
package_msg = ["PhysiologyAnalysis"]

function check_loaded_packages() 
     for package in package_msg
          println("$(package) is loaded")
     end
end
export check_loaded_packages

function __init__()
     @require RCall = "6f49c342-dc21-5d91-9882-a32aef131414" begin
          println("Loading R")
          export RCall
          push!(package_msg, "RCall")
     end

     @require FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341" begin
          include("Filtering/make_spectrum.jl")
          push!(package_msg, "FFTW")
     end

     #This is a good section to try using @Requires
     @require DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa" begin
          #using DifferentialEquations #For require, do we actually need to import this? 
          using DiffEqParamEstim, Optim
          include("Filtering/artifactRemoval.jl")
          export RCArtifact
          push!(package_msg, "DifferentialEquations")
     end
     #===============================Import all Datasheet tools==============================#
     #Only import if DataFrames has been loaded

     @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
          using RecipesBase
          include("Plotting/PhysRecipes.jl")
          push!(package_msg, "Plots(GR)")
     end

     @require ContinuousWavelets = "96eb917e-2868-4417-9cb6-27e7ff17528f" begin
          @require Wavelets = "29a6e085-ba6d-5f35-a997-948ac2efa89a" begin
               include("Filtering/wavelet_filtering.jl")
               export cwt_filter!, cwt_filter
               export dwt_filter!, dwt_filter
          end
     end
end

end