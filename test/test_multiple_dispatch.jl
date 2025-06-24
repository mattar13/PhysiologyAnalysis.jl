using Revise
using ElectroPhysiology
using PhysiologyAnalysis
using Statistics
import ElectroPhysiology: Experiment, TWO_PHOTON

# Test file for multiple dispatch functionality

println("Testing multiple dispatch for parameter loading functions...")

# This demonstrates how the multiple dispatch works:
# 
# 1. baseline_trace(data) - uses default parameters (existing function)
# 2. baseline_trace("parameters.json", data) - uses JSON parameters (new function)
#
# 1. process_rois(data) - uses default parameters (existing function)  
# 2. process_rois("parameters.json", data) - uses JSON parameters (new function)

println("\nMultiple dispatch examples:")
println("- baseline_trace(data) - uses default parameters")
println("- baseline_trace(\"parameters.json\", data) - uses JSON parameters")
println("- process_rois(data) - uses default parameters")
println("- process_rois(\"parameters.json\", data) - uses JSON parameters")

println("\nFunction signatures:")
println("- baseline_trace(trace::AbstractVector{T}; ...) - original function")
println("- baseline_trace(parameter_fn::String, data_img; ...) - parameterized version")
println("- process_rois(data::Experiment{TWO_PHOTON, T}; ...) - original function")
println("- process_rois(parameter_fn::String, data; ...) - parameterized version")

println("\nBenefits of multiple dispatch:")
println("- Clean API: same function name, different behavior based on arguments")
println("- Backward compatibility: existing code continues to work")
println("- Easy to use: just add parameter filename as first argument")
println("- Type safety: Julia automatically chooses the right method")

println("\nTest completed!") 