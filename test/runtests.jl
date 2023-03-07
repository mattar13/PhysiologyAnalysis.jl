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

    data_normalized = normalize(data)
    @test !isnothing(data_normalized)
end

@testset "Testing fitting functions" begin
    

end

@testset "Testing ERG analysis" begin
    data_filt = data_filter(data)
    resps = saturated_response(data_filt)
    @test maximum(resps) == -0.2988254587666856
end