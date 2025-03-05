const E_max = 30.0 # GeV
const m0_max = sqrt(2mp^2 + 2mp * E_max)

const model_zero = (;
    m0 = m0_max,
    masses_Pc = (4.312, 4.440, 4.457),
    widthes_Pc = (0.01, 0.02, 0.006),
    couplings_Pc = (0.1, 0.35, 0.08),
    scattlen_pp = 1.5,
    c_pp = 3.2 * cis(π / 3))

const ModelZero = typeof(model_zero)

k_pp(σ) = ThreeBodyDecays.breakup(sqrt(σ), mp, mp)
A_pp(σ; scattlen_pp) = 1 / (1 / scattlen_pp - 1im * k_pp(σ))

function ThreeBodyDecays.amplitude(model::ModelZero, σs)
    @unpack masses_Pc, widthes_Pc, couplings_Pc = model
    @unpack scattlen_pp, c_pp = model
    shapes = (BreitWigner(m, Γ) for (m, Γ) in zip(masses_Pc, widthes_Pc))
    σ1, σ2, σ3 = σs
    A = c_pp * A_pp(σ2; scattlen_pp) +
        sum(c * bw(σ1) for (c, bw) in zip(couplings_Pc, shapes)) +
        sum(c * bw(σ3) for (c, bw) in zip(couplings_Pc, shapes))
    return A
end

ThreeBodyDecays.masses(model::ModelZero) = ThreeBodyMasses(mp, mJψ, mp; model.m0)
