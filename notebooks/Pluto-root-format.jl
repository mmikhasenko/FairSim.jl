### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 8912962a-d8a2-11ef-1d62-3f9007747af6
begin
    cd(joinpath(@__DIR__, ".."))
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()

    using UnROOT
    using DataFrames
end

# ╔═╡ 69bd9219-d659-4790-ab69-19203f523fad
# Pkg.add("UnROOT")

# ╔═╡ a385c5c1-a13b-4626-9573-268015073dac
file_name = joinpath("data", "test_pluto_format.root");

# ╔═╡ 7d619620-116f-477d-acc6-412c76836f76
!isfile(file_name) && error("File is not found!")

# ╔═╡ e5a970c8-2d56-4980-91c7-c7fa0e488b3e
f = ROOTFile(file_name)

# ╔═╡ f6c68d27-7de3-4273-923b-aa020ac0c347
t = LazyTree(f, "data", ["Npart", "Impact", "Phi"]);

# ╔═╡ 8f05ac8e-d64c-446c-9867-affc7c6725b1
t

# ╔═╡ Cell order:
# ╠═8912962a-d8a2-11ef-1d62-3f9007747af6
# ╠═69bd9219-d659-4790-ab69-19203f523fad
# ╠═a385c5c1-a13b-4626-9573-268015073dac
# ╠═7d619620-116f-477d-acc6-412c76836f76
# ╠═e5a970c8-2d56-4980-91c7-c7fa0e488b3e
# ╠═f6c68d27-7de3-4273-923b-aa020ac0c347
# ╠═8f05ac8e-d64c-446c-9867-affc7c6725b1
