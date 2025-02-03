### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 0f61210a-a1db-11ef-2af3-e76254c194e4
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    # 
    using HadronicLineshapes
    using ThreeBodyDecays
    using LinearAlgebra
    using LaTeXStrings
    using DataFrames
    using Parameters
    using Setfield
    using QuadGK
    using PlutoUI
    using Plots
end

# ╔═╡ 07ed7deb-4829-4f67-bfa9-98479dcfe371
axis_labels = (;
	xlab = L"m^2(J/\psi\,p)\,\, [\textrm{GeV}^2]",
	ylab = L"m^2(J/\psi\,p)\,\, [\textrm{GeV}^2]");

# ╔═╡ dcdd4853-b8fa-4b8e-9547-508bda760339
theme(:boxed, size=(500,450))

# ╔═╡ 45789d1f-fdb4-41c5-b352-95f031355f5e
begin
	const mp = 0.984
	const mJψ = 3.09
end

# ╔═╡ c0eaacad-4be9-4f3f-8fe7-a6c08accd51d
const m0_max = sqrt(2mp*30)

# ╔═╡ 516dc354-0d3e-43f8-a613-bfe9226dacee
@bind m0 Slider(5.5:0.01:m0_max, default=m0_max, show_value=true)

# ╔═╡ 187f11a6-a8e5-4f09-8733-9d4799decf0f
ms = ThreeBodyMasses(mp, mJψ, mp; m0)

# ╔═╡ 206007e0-a462-4ec0-a090-320070e0efd8
σs_test = randomPoint(ms)

# ╔═╡ bdf603a8-d474-4388-bbdb-c42a16087d73
k_pp(σ) = ThreeBodyDecays.breakup(sqrt(σ), mp, mp)

# ╔═╡ d03292ab-9ec6-4434-9545-f06fdc4b215c
A_pp(σ; scattlen_pp) = 1/(1/scattlen_pp - 1im * k_pp(σ))

# ╔═╡ 5d7732c1-2ff3-4da6-b4ae-08d684571efa
function A(σs, pars)
	@unpack masses_Pc, widthes_Pc, couplings_Pc = pars
	@unpack scattlen_pp, c_pp = pars
	shapes = (BreitWigner(m, Γ) for (m,Γ) in zip(masses_Pc, widthes_Pc))
	σ1,σ2,σ3 = σs
	c_pp * A_pp(σ2; scattlen_pp) +
		sum(c*bw(σ1) for (c, bw) in zip(couplings_Pc, shapes)) +
		sum(c*bw(σ3) for (c, bw) in zip(couplings_Pc, shapes))
end

# ╔═╡ d15d368e-ad38-4b82-b41e-7174356842f4
const model = (; 
	masses_Pc = (4.312, 4.440, 4.457),
	widthes_Pc = (0.01, 0.02, 0.006),
	couplings_Pc = (0.1, 0.35, 0.08),
	scattlen_pp=1.5, # / GeV
	c_pp=3.2*cis(π/3))

# ╔═╡ 523445d2-a59e-41db-877c-15796494faff
begin
	plot(border31(ms), aspect_ratio=1)
	hline!(model.masses_Pc.^2 |> collect)
	vline!(model.masses_Pc.^2 |> collect)
	plot!(; axis_labels...)
	savefig("schematic_dalitz_ppJpsi.pdf")
	plot!()
end

# ╔═╡ 2020b904-59c1-4f4d-98a7-f0c8611ce348
begin
	plot(ms; iσx=1, iσy=3, aspect_ratio=1, grid_size=300,
		xlab="m²(J/ψ p) [GeV²]", ylab="m²(J/ψ p) [GeV²]") do σs
		a = A(σs, model)
		abs2(a)
	end
	plot!(; axis_labels...)
	savefig("dalitz_ppJpsi.pdf")
	plot!()
end

# ╔═╡ 25699e31-5214-4fd7-b463-601ebdfe5ffe
model_Pc = @set model.c_pp = 0.0im

# ╔═╡ d1692d47-67d5-4ec2-806a-c19fe14bb52d
model_pp = @set model.couplings_Pc = (0.0, 0.0, 0.0)

# ╔═╡ fe5a7dc6-e5b2-4a2d-90c3-4256aeb7aa2b
begin
	plot(leg=:bottom,
		ylab=L"\textrm{d}\Gamma / \textrm{d}m(J/\psi\,p)",
		size=(600,450))
	plot!(4.0, 6.9, fill=0, c=3, lc=:black, alpha=0.3, linealpha=1, lab="full") do e1
		integrand = projection_integrand(σs->abs2(A(σs, model)), ms, e1^2; k = 3)
		e1 * quadgk(integrand, 0, 1)[1]
	end
	plot!(4.0, 6.9, l=(:black, :dash), lw=0.4, lab="pp") do e1
		integrand = projection_integrand(σs->abs2(A(σs, model_pp)), ms, e1^2; k = 3)
		e1 * quadgk(integrand, 0, 1)[1]
	end
	lens!([4.25, 4.55], [0, 3200], inset=(1,bbox(0.0,0.0,0.5,0.5,:right, :top)))
	plot!(sp=2, yaxis=nothing, xlab=L"m(J/\psi\,p)\,\, [\textrm{GeV}]")
	savefig("projection_ppJpsi.pdf")
	plot!()
end

# ╔═╡ Cell order:
# ╠═0f61210a-a1db-11ef-2af3-e76254c194e4
# ╠═07ed7deb-4829-4f67-bfa9-98479dcfe371
# ╠═dcdd4853-b8fa-4b8e-9547-508bda760339
# ╠═c0eaacad-4be9-4f3f-8fe7-a6c08accd51d
# ╠═516dc354-0d3e-43f8-a613-bfe9226dacee
# ╠═45789d1f-fdb4-41c5-b352-95f031355f5e
# ╠═187f11a6-a8e5-4f09-8733-9d4799decf0f
# ╠═523445d2-a59e-41db-877c-15796494faff
# ╠═5d7732c1-2ff3-4da6-b4ae-08d684571efa
# ╠═206007e0-a462-4ec0-a090-320070e0efd8
# ╠═bdf603a8-d474-4388-bbdb-c42a16087d73
# ╠═d03292ab-9ec6-4434-9545-f06fdc4b215c
# ╠═d15d368e-ad38-4b82-b41e-7174356842f4
# ╠═2020b904-59c1-4f4d-98a7-f0c8611ce348
# ╠═25699e31-5214-4fd7-b463-601ebdfe5ffe
# ╠═d1692d47-67d5-4ec2-806a-c19fe14bb52d
# ╠═fe5a7dc6-e5b2-4a2d-90c3-4256aeb7aa2b
