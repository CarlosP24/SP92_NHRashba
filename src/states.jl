using JLD2
using Quantica
using Quantica: σ
σ0τz = σ(3, 0)
σzτ0 = σ(0, 3)
σyτz = σ(3, 2)
σ0τx = σ(1, 0)
σ0τ0 = σ(0, 0)
σyτ0 = σ(0, 2)
σxτ0 = σ(0, 1)
τz = σz = σ(3)
τx = σx = σ(1)
τy = σy = σ(2)
τ0 = σ0 = σ(0)
using ProgressMeter, Parameters
using LinearAlgebra, Arpack

include("builders/Wire.jl")
include("builders/leads.jl")

include("models/systems.jl")

include("calculations/States.jl")

##
name = ARGS[1]
Vz = parse(Float64, ARGS[2])
nev = parse(Int, ARGS[3])
    
res = calc_states(name, Vz, nev)
save(res.path, "res", res)

