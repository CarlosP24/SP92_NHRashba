using Pkg
Pkg.activate("plots")
Pkg.instantiate()
using CairoMakie, JLD2, Parameters, Quantica, Interpolations
include("../src/builders/Wire.jl")
include("plotters.jl")