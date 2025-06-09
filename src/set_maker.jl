using Quantica, Parameters
include("builders/Kitaev.jl")
include("builders/Wire.jl")
include("builders/Filter.jl")
include("builders/leads.jl")

include("models/systems.jl")
include("models/Filters.jl")

keyword = ARGS[1]

open("sets/$(keyword)", "w") do io
    [println(io, m) for m in filter(s -> startswith(s, keyword), keys(systems))]
end