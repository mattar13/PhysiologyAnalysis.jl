module PhysiologyAnalysis

using Requires
# The top level is the ElectroPhysiology package. These are not imported into the workspace
using Dates

using ElectroPhysiology
import ElectroPhysiology: Experiment, readABF, parseABF
import ElectroPhysiology: now, year, month, day, hour, minute, second
import ElectroPhysiology: TWO_PHOTON
import ElectroPhysiology: readABFInfo
#= Packages used for fitting data ====================================#
using LsqFit #Used for fitting amplification, Intensity Response, and Resistance Capacitance models

#= Packages used for Analyzing data ==================================#
import Polynomials as PN #used for fitting and stats

using StatsBase #Used for mean, std, and SEM functions.
using Statistics
using HypothesisTests
import Statistics.mean

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
include("Analysis/ERGAnalysis/ERGAnalysis.jl")
#export calculate_basic_stats
export saturated_response, dim_response
export minima_to_peak, time_to_peak
export percent_recovery_interval #This finds the dominant time constant
export recovery_time_constant #This finds the recovery time constant
export integral #This finds the integration time
export get_response
export calculate_threshold

#using JLD2 #Maybe this should be added to Requires.jl
include("Analysis/WaveAnalysis/thresholding.jl")
export calculate_threshold

include("Analysis/WaveAnalysis/TimescaleAnalysis.jl")
export get_timestamps, extract_interval
export max_interval_algorithim, timeseries_analysis

# These functions are used by the base
#This file contains things like extraction and convienance functions
function set_calibration_path(pathname::String ;path = "$(homepath)/Datasheets/calibration.txt")
     open(path, "w") do file
          write(file, pathname)
     end
end

homepath = joinpath(splitpath(pathof(PhysiologyAnalysis))[1:end-1]...)

calibration_path() = read("$(homepath)/Datasheets/calibration.txt", String)

#3) Import all Datasheet tools ===========================================================#
function __init__()
     @require OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed" begin
          using .OrdinaryDiffEq
          @require ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78" begin
               println("Loading differential equations and nose fitting")
               using .ModelingToolkit
               include("Fitting/NoseModel.jl")
               export findNosePeak
          end
     end

     #Eventually we should massively restructure this
     @require XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0" begin
          #println(XLSX exported)
          println("Loading dataframes and file opening")
          using .XLSX
          @require DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0" begin
               using .DataFrames
               @require Query = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1" begin
                    using .Query 

                    include("Datasheets/FilePathExtraction.jl")
                    export traverse_root
                    export getABF_datetime
               end
          end
     end

     @require FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549" begin
          using .FileIO
          @require Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0" begin
               using .Images
               @require ImageView = "86fae568-95e7-573e-a6b2-d8a6b900c9ef" begin
                    using .ImageView
                    include("Analysis/ImagingAnalysis/PixelExtraction.jl")
                    println("Imported necessary things")
                    export zProject, frameAverage 
                    export normalize, binarize
                    export findROIcentroid
                    @require ImageFiltering = "6a3955dd-da59-5b1f-98d4-e7296123deb5" begin
                         include("Analysis/ImagingAnalysis/deltaF.jl")
                         export deltaF, deltaF_F, roll_mean
                    end
              end
          end
     end

     @require Flux = "587475ba-b771-5e3f-ad9e-33799f191a9c" begin
          using .Flux
          #include("Analysis/ImagingAnalysis/CellPose_port.jl")
          #export build_model
     end

     #We want to add something for PyCall 
end

end