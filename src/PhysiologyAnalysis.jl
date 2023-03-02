module PhysiologyAnalysis

# The top level is the ElectroPhysiology package
using ElectroPhysiology
import ElectroPhysiology: Experiment, readABF, parseABF
using DSP #Used for lowpass, highpass, EI, and notch filtering
import Polynomials as PN #

#export some basic functions from 
export readABF, parseABF
#This package does 4 things: 

#1)Filter ====================================================================================#
include("Filtering/filtering.jl")
export filter_data, filter_data!
export rolling_mean
export normalize, normalize!

include("Filtering/filteringPipelines.jl")
export data_filter!, data_filter

#2) Fitting ============================================================================#
using LsqFit #Used for fitting amplification, Intensity Response, and Resistance Capacitance models
using Distributions
include("Fitting/Models.jl")
export HILL_MODEL
export amplification
export curve_fit #curve fitting from LsqFit
export IR_curve

#3) Data anlysis ========================================================================#
using Statistics, StatsBase #These functions use R functions as well as StatsBase
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

#4) Import all Datasheet tools ===========================================================#
#Only import if DataFrames has been loaded
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

#5) Plotting utilities will be loaded in automatically ==============================================#
import PyPlot
import PyPlot.plt #All the base utilities for plotting
export plt
import PyPlot.matplotlib
import PyCall as py #This allows us to use Python to call somethings 
import PyCall: @pyimport, PyObject

include("Plotting/DefaultSettings.jl") #This requires PyPlot
include("Plotting/PlottingUtilities.jl")
include("Plotting/PhysPyPlot.jl")
export plot_experiment, plot_experiment_fit

include("Plotting/DatasheetPlotting.jl")
export plot_ir_scatter, plot_ir_fit, plot_IR


#include("Filtering/make_spectrum.jl")
#using ContinuousWavelets, Wavelets
#include("Filtering/wavelet_filtering.jl")
#export cwt_filter!, cwt_filter
#export dwt_filter!, dwt_filter

end