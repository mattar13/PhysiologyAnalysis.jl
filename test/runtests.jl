using ElectroPhysiology
using PhysiologyAnalysis
using Test

test_file = raw"to_analyze.abf"
data = readABF(test_file) |> data_filter

@testset "Testing ElectroPhysiology" begin
    @test !isnothing(data)
    @test isa(data, ElectroPhysiology.Experiment)
end

@testset "Testing ERG analysis" begin
    resps = saturated_response(data)
    #println(maximum(resps))
    #@test maximum(resps) == -0.2988254587666856
end