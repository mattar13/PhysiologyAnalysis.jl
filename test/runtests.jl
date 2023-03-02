using ElectroPhysiology
using PhysiologyAnalysis
using Test

test_file = raw"to_analyze.abf"
data = readABF(test_file)
@testset "Testing ElectroPhysiology" begin
    @test !isnothing(data)
    @test isa(data, ElectroPhysiology.Experiment)
end

@testset "Testing data filtering" begin
    data_filtered = filter_data(data)
    @test !isnothing(data_filtered)

    
end
