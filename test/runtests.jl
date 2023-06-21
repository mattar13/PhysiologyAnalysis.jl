using ElectroPhysiology
using PhysiologyAnalysis
using Test

test_file = raw"to_analyze.abf"
data = readABF(test_file) |> data_filter

@testset "Testing ERG analysis" begin
    include("testAnalysis.jl")
end