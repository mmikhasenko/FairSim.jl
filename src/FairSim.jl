module FairSim

using HadronicLineshapes
using ThreeBodyDecays
using LinearAlgebra
using DataFrames
using Parameters
using Setfield
using QuadGK


include("constants.jl")

export masses
export amplitude
include("model_zero.jl")

end # module FairSim
