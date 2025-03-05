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
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    # 
    using HadronicLineshapes
    using ThreeBodyDecays
    using LinearAlgebra
    using LaTeXStrings
    using Measurements
    using DataFrames
    using Parameters
    using Statistics
    using Setfield
    using QuadGK
    using PlutoUI
    using Plots
	# 
	using FairSim
	import FairSim: m0_max
end

# ╔═╡ 783954c5-0123-4c2b-a28d-1aeb7eb19218
const plot_folder = joinpath(@__DIR__, "..", "plots")

# ╔═╡ 07ed7deb-4829-4f67-bfa9-98479dcfe371
axis_labels = (;
	xlab = L"m^2(J/\psi\,p)\,\, [\textrm{GeV}^2]",
	ylab = L"m^2(J/\psi\,p)\,\, [\textrm{GeV}^2]");

# ╔═╡ dcdd4853-b8fa-4b8e-9547-508bda760339
theme(:boxed, size=(500,450))

# ╔═╡ 516dc354-0d3e-43f8-a613-bfe9226dacee
@bind m0 Slider(5.5:0.01:m0_max, default=m0_max, show_value=true)

# ╔═╡ 45789d1f-fdb4-41c5-b352-95f031355f5e
const model = (;
    m0,
    masses_Pc = (4.312, 4.440, 4.457),
    widthes_Pc = (0.01, 0.02, 0.006),
    couplings_Pc = (0.1, 0.35, 0.08),
    scattlen_pp = 1.5,
    c_pp = 3.2 * cis(π / 3))

# ╔═╡ 187f11a6-a8e5-4f09-8733-9d4799decf0f
ms = masses(model)

# ╔═╡ 523445d2-a59e-41db-877c-15796494faff
begin
	plot(border31(ms), aspect_ratio=1)
	hline!(model.masses_Pc.^2 |> collect)
	vline!(model.masses_Pc.^2 |> collect)
	plot!(; axis_labels...)
	savefig(joinpath(plot_folder, "schematic_dalitz_ppJpsi.pdf"))
	plot!()
end

# ╔═╡ 2020b904-59c1-4f4d-98a7-f0c8611ce348
begin
	plot(ms; iσx=1, iσy=3, aspect_ratio=1, grid_size=300,
		xlab="m²(J/ψ p) [GeV²]", ylab="m²(J/ψ p) [GeV²]") do σs
		a = amplitude(model, σs)
		abs2(a)
	end
	plot!(; axis_labels...)
	savefig(joinpath(plot_folder, "dalitz_ppJpsi.pdf"))
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
		integrand = projection_integrand(σs->abs2(amplitude(model, σs)), ms, e1^2; k = 3)
		e1 * quadgk(integrand, 0, 1)[1]
	end
	plot!(4.0, 6.9, l=(:black, :dash), lw=0.4, lab="pp") do e1
		integrand = projection_integrand(σs->abs2(amplitude(model_pp, σs)), ms, e1^2; k = 3)
		e1 * quadgk(integrand, 0, 1)[1]
	end
	lens!([4.25, 4.55], [0, 3200], inset=(1,bbox(0.0,0.0,0.5,0.5,:right, :top)))
	plot!(sp=2, yaxis=nothing, xlab=L"m(J/\psi\,p)\,\, [\textrm{GeV}]")
	savefig(joinpath(plot_folder, "projection_ppJpsi.pdf"))
	plot!()
end

# ╔═╡ 1246dd11-8168-4c04-97f9-eb6495edbd76
md"""
## Fit fractions
"""

# ╔═╡ 68a4140d-1288-4a59-ae93-5bbec9757a9e
const model29 = let
	E_max = 29.0
	mp = FairSim.mp
	model_zero = FairSim.model_zero
	@set model_zero.m0 = sqrt(2mp^2 + 2mp * E_max)
end;

# ╔═╡ 0b553929-7661-427d-a3d9-84268e509e48
const model29_Pc = @set model29.c_pp = 0.0im;

# ╔═╡ 3dc2d704-cc12-45cb-9007-5ddba15d4bfd
const model29_pp = @set model29.couplings_Pc = (0.0, 0.0, 0.0);

# ╔═╡ 1f0511f1-6e78-4630-88c6-ed3e90449490
function numerical_integral(model, phsp)
	_i = map(σs->abs2(amplitude(model, σs)), phsp)
	mean(_i) ± std(_i) / sqrt(length(_i))
end

# ╔═╡ 2ee4b4c3-6471-4856-98b9-7c29c5d141de
integral_computations = let
	nDraft = 100_000
	ms = masses(model29)
	_data = mapslices(rand(nDraft,2); dims=2) do y
		y2σs(y, ms; k=1)
	end[:,1];
	filter!(_data) do σs
		isphysical(σs, ms)
	end
	σs0 = randomPoint(ms)
	I_all = numerical_integral(model29, _data)
	I_pp = numerical_integral(model29_pp, _data)
	I_Pcc = map(1:3) do i
		c = model29.couplings_Pc
		_model_pure = @set model29_Pc.couplings_Pc = ntuple(x->0.0, length(c))
		_model = @set _model_pure.couplings_Pc[i] = c[i]
		numerical_integral(_model, _data)
	end
	# 
	(; I_all, I_pp, I_Pcc)
end

# ╔═╡ ab049f0d-d483-44f2-bb59-d4c3c63c35c5
fit_fractions = (
	:pp => integral_computations.I_pp / integral_computations.I_all .* 100,
	([:Pcc4312, :Pcc4440, :Pcc4457] .=> (integral_computations.I_Pcc ./ integral_computations.I_all .* 100))...
)

# ╔═╡ Cell order:
# ╠═0f61210a-a1db-11ef-2af3-e76254c194e4
# ╠═783954c5-0123-4c2b-a28d-1aeb7eb19218
# ╠═07ed7deb-4829-4f67-bfa9-98479dcfe371
# ╠═dcdd4853-b8fa-4b8e-9547-508bda760339
# ╠═516dc354-0d3e-43f8-a613-bfe9226dacee
# ╠═45789d1f-fdb4-41c5-b352-95f031355f5e
# ╠═187f11a6-a8e5-4f09-8733-9d4799decf0f
# ╠═523445d2-a59e-41db-877c-15796494faff
# ╠═2020b904-59c1-4f4d-98a7-f0c8611ce348
# ╠═25699e31-5214-4fd7-b463-601ebdfe5ffe
# ╠═d1692d47-67d5-4ec2-806a-c19fe14bb52d
# ╟─1246dd11-8168-4c04-97f9-eb6495edbd76
# ╠═68a4140d-1288-4a59-ae93-5bbec9757a9e
# ╠═0b553929-7661-427d-a3d9-84268e509e48
# ╠═3dc2d704-cc12-45cb-9007-5ddba15d4bfd
# ╠═1f0511f1-6e78-4630-88c6-ed3e90449490
# ╠═2ee4b4c3-6471-4856-98b9-7c29c5d141de
# ╟─ab049f0d-d483-44f2-bb59-d4c3c63c35c5
# ╠═fe5a7dc6-e5b2-4a2d-90c3-4256aeb7aa2b
