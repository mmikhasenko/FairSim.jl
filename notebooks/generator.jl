### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ bc2de5d4-b3d2-11ef-2944-0986d8aaa462
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    # 
    using HadronicLineshapes: BreitWigner
    using ThreeBodyDecays
    using LinearAlgebra
    using LaTeXStrings
    using FourVectors
	using DataFrames
    using Parameters
    using Setfield
    using QuadGK
    using PlutoUI
    using Plots
	using CSV
	# 
	using FairSim
end

# ╔═╡ 5b12360b-5b5a-4163-b0ec-6018d67de0da
md"""
# Importance sampling

In this notebook a four-vectors are created from set of Mandelstam invariant variables. The vectors are boosted to the lab frame using **fixed target** setup with $30\,$GeV proton beam on a proton target.

We use the `FairSim.model_zero` to weight the events. It generates a matrix element for $J/\psi$ production

$p p \to p p J/\psi$

using a simple BW model.
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

# ╔═╡ 1d034ecb-afb8-40a3-b902-f25f6d18475d
const model = FairSim.model_zero

# ╔═╡ a9cd25ec-17d4-4cd3-be2d-ab0cb1ab87dc
const ms = masses(model)

# ╔═╡ 6c4c35fe-77c9-4719-8bfd-da7580061e65
function propagate_to_lab(p, euler_angles, γ_lab = FairSim.E_max/ms.m0)
	@unpack α,β,γ = euler_angles
	p_rotated = p |> Rz(γ) |> Ry(β) |> Rz(α) 
	p_lab = p_rotated |> Bz(γ_lab)
	return p_lab
end

# ╔═╡ c0c594bf-f92d-4d1e-b708-37c48a24d68b
function decay_a_particle(p; decay_angles::NamedTuple{(:cosθ, :ϕ)}, m1, m2)
	# parameters of the mother moving frame
	γ = boost_gamma(p)
	m = mass(p)
	@unpack cosθ, ϕ = spherical_coordinates(p)

	# get vectors in rest frame
	k = breakup(m,m1,m2)
	sinθ_decay = sqrt(1-decay_angles.cosθ^2)
	k_vec = (
		k * sinθ_decay * cos(decay_angles.ϕ),
		k * sinθ_decay * sin(decay_angles.ϕ), k * decay_angles.cosθ)
	# 
	# get four-vectors in lab frame
	p1 = FourVector(k_vec...; M = m1)
	p2 = FourVector((0 .- k_vec)...; M = m2)
	(p1, p2) .|> Bz(γ) .|> Ry(acos(cosθ)) .|> Rz(ϕ)
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
	propagate_to_lab.(pv_0, (; α=0.1, β=0.2, γ=0.3) |> Ref)
end;

# ╔═╡ 0f4d1022-c1e3-424f-8df8-e59bb7a01477
let # test
	p = four_vectors_test[1]
	p1, p2 = decay_a_particle(p; decay_angles=(cosθ=0.6, ϕ=0.2), m1=0.2, m2=0.2)
	@assert p ≈ (p1 + p2)
end

# ╔═╡ 97074042-4620-43e1-af12-0b7d189226c2
@assert sum(four_vectors_test)[4] ≈ FairSim.E_max

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
function kinematic_mapping(y; k)
	_σs = x2σs(y[1:2], ms; k)
	pv_0 = reference_four_vectors(_σs, ms)	
	# 
	euler_angles = (; 
		α=π*(2y[3]-1),
		β=acos(2y[4]-1),
		γ=π*(2y[5]-1)
		)
	#
	propagate_to_lab.(pv_0, Ref(euler_angles))
end

# ╔═╡ b374abc9-1934-43a8-a113-3b8443ae67b5
to_vector(p) = (@unpack px,py,pz,E = p; [px,py,pz,E])

# ╔═╡ 5acf4dfc-85d2-450a-9c0d-03fa73bec705
const nMC = 100_000;

# ╔═╡ e26db898-6dcc-411b-a66c-d98a6ec62150
const k_preference = 3

# ╔═╡ cfad287a-7303-43e0-ba05-d96407e6580b
data_x = let
	sample0 = rand(nMC,5);
	mapslices(sample0; dims=2) do y
		pv = kinematic_mapping(y; k=k_preference)
		NamedTuple{(:p1_xyzE,:p2_xyzE,:p3_xyzE)}(pv)
	end[:,1]
end |> DataFrame;

# ╔═╡ 8b532d54-394b-496e-9838-79cca12e65af
data_y = transform(data_x, [:p1_xyzE, :p2_xyzE, :p3_xyzE] => ByRow() do p1, p2, p3
	σ3 = mass2(p1+p2)
	σ1 = mass2(p2+p3)
	σs = Invariants(ms; σ1, σ3)
	a = amplitude(model, σs)
	# 
	k = k_preference
	σk = σs[k]
	weight_phsp = breakup_ij(σk, ms; k)*breakup_Rk(σk, ms; k) / sqrt(σk)*ms.m0;
	weight = abs2(a) * weight_phsp
	# 
	(; weight, σs)
end => AsTable);

# ╔═╡ 1297b8ea-41f8-4e3c-836d-276fe6ace8ef
md"""
### Importance Sampling
"""

# ╔═╡ 5cfd9476-ed34-4f54-ba2f-2ed3f9b49a39
function importance_sampling(weights)
	max_w = maximum(weights)
	n = length(weights)
	if_to_take = weights .> rand(n) .* max_w
	return if_to_take
end

# ╔═╡ e6e92ca6-d368-48ec-abae-2bb32644daec
data_sampled = data_y[importance_sampling(data_y.weight), :];

# ╔═╡ 391df748-973c-4174-85aa-8bf557d785e6
let bins=100
	plot(size=(700,300), grid=(1,2),
		histogram2d(getindex.(data_y.σs, 1), getindex.(data_y.σs, 3); bins),
		histogram2d(getindex.(data_sampled.σs, 1), getindex.(data_sampled.σs, 3); bins))
end

# ╔═╡ bf0eed65-4848-4315-a1c4-73210f2b3cf5
data_with_muons_sampled = transform(data_sampled, :p2_xyzE => ByRow() do p
	pmup_xyzE, pmum_xyzE = decay_a_particle(p;
		decay_angles=(cosθ=0.6, ϕ=0.2), m1=FairSim.mμ, m2=FairSim.mμ)
	(; pmup_xyzE, pmum_xyzE)
end => AsTable);

# ╔═╡ f6927ce9-6c1e-413d-bece-78c813d801c5
CSV.write(joinpath(@__DIR__, "..", "data", "momenta.csv"), 
	select(data_with_muons_sampled, [:p1_xyzE, :pmup_xyzE, :pmum_xyzE, :p3_xyzE]))

# ╔═╡ 8a3240a7-453e-489e-8eab-6194de710c52
md"""
### sample size: $(size(data_sampled, 1))
"""

# ╔═╡ Cell order:
# ╟─5b12360b-5b5a-4163-b0ec-6018d67de0da
# ╠═bc2de5d4-b3d2-11ef-2944-0986d8aaa462
# ╠═089ef92b-4c64-44ef-9d44-77f8ae0aba14
# ╠═1d034ecb-afb8-40a3-b902-f25f6d18475d
# ╠═a9cd25ec-17d4-4cd3-be2d-ab0cb1ab87dc
# ╠═6c4c35fe-77c9-4719-8bfd-da7580061e65
# ╠═c0c594bf-f92d-4d1e-b708-37c48a24d68b
# ╠═0f4d1022-c1e3-424f-8df8-e59bb7a01477
# ╟─c58ddd76-08cc-40c6-8c8b-c4f025b87bf9
# ╠═7941cb28-d2c9-4a4c-9db1-565c755fbc91
# ╠═dcc69136-2ff2-46d8-ab96-c2888f0610d0
# ╠═97074042-4620-43e1-af12-0b7d189226c2
# ╠═8b611485-2b76-4238-a689-48e2689181cd
# ╟─99801d70-a388-4df2-a2b8-e5510200bc00
# ╠═75d147ba-da11-4c0f-9736-03cd3d86c45c
# ╠═b374abc9-1934-43a8-a113-3b8443ae67b5
# ╠═5acf4dfc-85d2-450a-9c0d-03fa73bec705
# ╠═e26db898-6dcc-411b-a66c-d98a6ec62150
# ╠═cfad287a-7303-43e0-ba05-d96407e6580b
# ╠═8b532d54-394b-496e-9838-79cca12e65af
# ╟─1297b8ea-41f8-4e3c-836d-276fe6ace8ef
# ╠═5cfd9476-ed34-4f54-ba2f-2ed3f9b49a39
# ╠═e6e92ca6-d368-48ec-abae-2bb32644daec
# ╠═391df748-973c-4174-85aa-8bf557d785e6
# ╠═bf0eed65-4848-4315-a1c4-73210f2b3cf5
# ╠═f6927ce9-6c1e-413d-bece-78c813d801c5
# ╟─8a3240a7-453e-489e-8eab-6194de710c52
