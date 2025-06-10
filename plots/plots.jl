using Pkg
Pkg.activate("plots")
Pkg.instantiate()
using CairoMakie, JLD2, Parameters, Quantica
include("../src/builders/Kitaev.jl")
include("../src/builders/Wire.jl")
include("../src/builders/Filter.jl")
include("plotters.jl")




