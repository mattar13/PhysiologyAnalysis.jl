using Plots, RecipesBase 
using Colors, StatsPlots
import PyPlot as plt #All the base utilities for plotting
import PyPlot.matplotlib

include("DefaultSettings.jl") #This imports all the default plotting settings for PyPlot
include("PlottingUtilities.jl") #This imports all the plotting utilites
include("PhysRecipes.jl")
#include("PhysPyPlot.jl") #Imports all the PyPlot utilities

export plt #Export plotting utilities
export plot_experiment