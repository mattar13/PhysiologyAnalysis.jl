module PhysiologyAnalysis

using Requires
# The top level is the ElectroPhysiology package. These are not imported into the workspace
using Dates

using ElectroPhysiology
import ElectroPhysiology: Experiment, readABF, parseABF
import ElectroPhysiology: eachtrial, eachchannel
import ElectroPhysiology: now, year, month, day, hour, minute, second
import ElectroPhysiology: TWO_PHOTON
import ElectroPhysiology: readABFInfo, getABF_datetime

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
using LsqFit
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

include("Analysis/WaveAnalysis/thresholding.jl")
export calculate_threshold

using Peaks #Use this as a better way to find peaks
include("Analysis/WaveAnalysis/TimescaleAnalysis.jl")
export get_timestamps, extract_interval
export max_interval_algorithim, timeseries_analysis
export findmaxima
export crosscor, crosscor!

include("Analysis/WholeCellAnalysis/passive_analysis.jl")
export calculate_baseline, calculate_peak 
export calculate_resistance, calculate_capacitance
export extract_timepoint

# These functions are used by the base
#This file contains things like extraction and convienance functions
function set_calibration_path(pathname::String ;path = "$(homepath)/Datasheets/calibration.txt")
     open(path, "w") do file
          write(file, pathname)
     end
end

using FileIO, Images
include("Analysis/ImagingAnalysis/PixelExtraction.jl")
export zProject, frameAverage 
export normalize, binarize
export findROIcentroid

using SparseArrays, OffsetArrays, ImageFiltering
include("Analysis/ImagingAnalysis/DeltaFF.jl")
export baseline_trace, baseline_stack

include("Analysis/ImagingAnalysis/ParametricFit.jl")
export single_stim_model

include("Analysis/ImagingAnalysis/ROIAnalysis.jl")
export findROIcentroids
export pixel_splits

using Interpolations
include("Analysis/Stats.jl")
export cor_xy

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
                    export parse_cell_details

                    include("Datasheets/DataSheetCreation.jl")
                    export create2PDataSheet
                    export save2PDataSheet, open2PDataSheet

                    include("Datasheets/DataSheetModify.jl")
                    export expand_dates

                    include("Datasheets/DataSheetAnalysis.jl")
                    export pair_experiments!, IV_analysis!
               end
          end
     end

     @require PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin       
          using .PyCall
          @require Conda = "8f4d0f93-b110-5947-807f-2305c1781a2d" begin
               using .Conda
               include("Analysis/ImagingAnalysis/CellPose.jl")
               println("CellPose port loaded")
               export cellpose_model
          end
     end

     #We want to add something for PyCall 
end

end