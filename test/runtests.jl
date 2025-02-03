using ThreeBodyDecays
using Parameters
using FairSim
using Random
using Test

Random.seed!(1234)


model = FairSim.model_zero

@testset "Amplitude Tests" begin
    ms = masses(model)
    σs = randomPoint(ms)
    amp = amplitude(model, σs)
    @show amp
    @test amp isa Complex
    @test amp ≈ -1.2517217630985085 + 1.4906520871784779im
end
