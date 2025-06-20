using Quantica, Parameters
include("builders/Wire.jl")
include("builders/leads.jl")

include("models/systems.jl")
include("parallelizers/pgeneric.jl")
include("parallelizers/pspectrum.jl")

include("calculations/Spectrum.jl")
include("calculations/Conductance.jl")

keyword = ARGS[1]

open("sets/$(keyword)", "w") do io
    [println(io, m) for m in filter(s -> startswith(s, keyword), keys(systems))]
end