using ThreeBodyDecays
using FairSim
using Test

@testset "Amplitude Tests" begin
    model = FairSim.model_zero
    ms = masses(model)
    # Fix the physical point rather than tying the regression to an RNG stream.
    σs = Invariants(ms; σ1 = 25.0, σ2 = 15.0)
    @test isphysical(σs, ms)
    amp = amplitude(model, σs)
    @test amp isa Complex

    # Independent closed form: equal-mass pp momentum and constant-width
    # Breit-Wigner propagators (the Pc shapes have massless daughters and L=0).
    k = sqrt(σs.σ2 / 4 - FairSim.mp^2)
    expected = model.c_pp / (inv(model.scattlen_pp) - im * k)
    for (m, Γ, c) in zip(model.masses_Pc, model.widthes_Pc, model.couplings_Pc)
        expected += c / (m^2 - σs.σ1 - im * m * Γ)
        expected += c / (m^2 - σs.σ3 - im * m * Γ)
    end
    @test amp ≈ expected

    # The two proton-pair assignments enter with identical Pc couplings.
    swapped = Invariants(ms; σ1 = σs.σ3, σ2 = σs.σ2)
    @test amplitude(model, swapped) ≈ amp
end
