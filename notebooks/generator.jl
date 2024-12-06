### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ bc2de5d4-b3d2-11ef-2944-0986d8aaa462
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec(url = "https://github.com/JuliaHEP/LorentzVectorBase.jl"),
		Pkg.PackageSpec(url = "https://github.com/mmikhasenko/FourVectors.jl"),
		Pkg.PackageSpec("ThreeBodyDecays"),
		Pkg.PackageSpec("Parameters"),
		Pkg.PackageSpec("DataFrames"),
		Pkg.PackageSpec("CSV")
	])

	using Parameters
	using DataFrames
	using FourVectors
	using Plots
	using ThreeBodyDecays
	using CSV
end

# ╔═╡ 5b12360b-5b5a-4163-b0ec-6018d67de0da
md"""
# Four-vectors from rest fram to lab frame

In this notebook a four-vectors are created from set of Mandelstam invariant variables. The vectors are boosted to the lab frame using **fixed target** setup with $30\,$GeV proton beam on a proton target.

$p p \to p p J/\psi$
"""

# ╔═╡ 089ef92b-4c64-44ef-9d44-77f8ae0aba14
function reference_four_vectors(σs, ms)
	m1sq,m2sq,m3sq,s = ms^2
	@unpack m0 = ms
	@unpack σ1,σ2,σ3 = σs
	# 
	E1 = (s+m1sq-σ1) / (2m0)
	E2 = (s+m2sq-σ2) / (2m0)
	E3 = (s+m3sq-σ3) / (2m0)
	# 
	p1 = sqrt(Kallen(s,m1sq,σ1)) / (2m0)
	p2 = sqrt(Kallen(s,m2sq,σ2)) / (2m0)
	p3 = sqrt(Kallen(s,m3sq,σ3)) / (2m0)
	# 
	cosθ = -cosζ(wr(1,2,0), σs, ms^2) # pi-angle(1,2) in rest of 0
	sinθ = sqrt(1-cosθ^2)
	# 
	FourVector(0,0,-p1, E1),
	FourVector(p2*sinθ,0,p2*cosθ, E2),
	FourVector(-p2*sinθ,0,p1-p2*cosθ, E3)
end

# ╔═╡ 2c092a3d-4020-41fe-8522-a72a1b9f1e7a
begin
	const mp = 0.938; # GeV
	const mψ = 3.09; # GeV
end

# ╔═╡ b9419922-0dc7-4e48-9591-74d6a5ef8a57
const E_beam = 30; # GeV

# ╔═╡ e1500a9b-47f7-46cd-9163-9d6e814823d5
sqrt(2mp^2+2E_beam*mp)

# ╔═╡ a9cd25ec-17d4-4cd3-be2d-ab0cb1ab87dc
ms = ThreeBodyMasses(mp, mψ, mp; m0 = sqrt(2mp^2+2E_beam*mp))

# ╔═╡ 6c4c35fe-77c9-4719-8bfd-da7580061e65
function propagate_to_lab(pv_0, euler_angles, γ_lab = E_beam/ms.m0)
	@unpack α,β,γ = euler_angles
	pv_rf = pv_0 .|> Rz(γ) .|> Ry(β) .|> Rz(α) 
	pv_rf .|> Bz(γ_lab)
end

# ╔═╡ c58ddd76-08cc-40c6-8c8b-c4f025b87bf9
md"""
## Tests of implementation
"""

# ╔═╡ 7941cb28-d2c9-4a4c-9db1-565c755fbc91
σs_test = randomPoint(ms)

# ╔═╡ dcc69136-2ff2-46d8-ab96-c2888f0610d0
four_vectors_test = let
	pv_0 = reference_four_vectors(σs_test, ms)	
	propagate_to_lab(pv_0, (; α=0.1, β=0.2, γ=0.3))
end;

# ╔═╡ 97074042-4620-43e1-af12-0b7d189226c2
@assert sum(four_vectors_test)[4] ≈ E_beam

# ╔═╡ 8b611485-2b76-4238-a689-48e2689181cd
@assert all(
	(mass2(sum(four_vectors_test[[2,3]])),
	mass2(sum(four_vectors_test[[3,1]])),
	mass2(sum(four_vectors_test[[1,2]]))) .≈ Tuple(collect(σs_test))
	)

# ╔═╡ 99801d70-a388-4df2-a2b8-e5510200bc00
md"""
## Generate data
"""

# ╔═╡ 75d147ba-da11-4c0f-9736-03cd3d86c45c
function kinematic_mapping(y)
	_σs = x2σs(y[1:2], ms; k=1)
	pv_0 = reference_four_vectors(_σs, ms)	
	# 
	euler_angles = (; 
		α=π*(2y[3]-1),
		β=acos(2y[4]-1),
		γ=π*(2y[5]-1)
		)
	#
	propagate_to_lab(pv_0, euler_angles)
end

# ╔═╡ 799bacec-5fad-4e21-876a-645a401c285f
kinematic_mapping(rand(5))

# ╔═╡ b374abc9-1934-43a8-a113-3b8443ae67b5
to_vector(p) = (@unpack px,py,pz,E = p; [px,py,pz,E])

# ╔═╡ 8a3240a7-453e-489e-8eab-6194de710c52
begin # write to disc
	data = map(1:10) do _
		pv = kinematic_mapping(rand(5))
		NamedTuple{(:p1_xyzE,:p2_xyzE,:p3_xyzE)}(pv)
	end |> DataFrame
	CSV.write(joinpath(@__DIR__, "..", "data", "momenta.csv"), data)
end;

# ╔═╡ Cell order:
# ╟─5b12360b-5b5a-4163-b0ec-6018d67de0da
# ╠═bc2de5d4-b3d2-11ef-2944-0986d8aaa462
# ╠═089ef92b-4c64-44ef-9d44-77f8ae0aba14
# ╠═2c092a3d-4020-41fe-8522-a72a1b9f1e7a
# ╠═b9419922-0dc7-4e48-9591-74d6a5ef8a57
# ╠═e1500a9b-47f7-46cd-9163-9d6e814823d5
# ╠═a9cd25ec-17d4-4cd3-be2d-ab0cb1ab87dc
# ╠═6c4c35fe-77c9-4719-8bfd-da7580061e65
# ╟─c58ddd76-08cc-40c6-8c8b-c4f025b87bf9
# ╠═7941cb28-d2c9-4a4c-9db1-565c755fbc91
# ╠═dcc69136-2ff2-46d8-ab96-c2888f0610d0
# ╠═97074042-4620-43e1-af12-0b7d189226c2
# ╠═8b611485-2b76-4238-a689-48e2689181cd
# ╟─99801d70-a388-4df2-a2b8-e5510200bc00
# ╠═75d147ba-da11-4c0f-9736-03cd3d86c45c
# ╠═799bacec-5fad-4e21-876a-645a401c285f
# ╠═b374abc9-1934-43a8-a113-3b8443ae67b5
# ╠═8a3240a7-453e-489e-8eab-6194de710c52
