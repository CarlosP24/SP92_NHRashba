using JLD2
@everywhere begin
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

    include("builders/Kitaev.jl")
    include("builders/Wire.jl")
    include("builders/Filter.jl")
    include("builders/leads.jl")

    include("models/systems.jl")
    include("parallelizers/pgeneric.jl")
    include("parallelizers/pspectrum.jl")

    include("calculations/Spectrum.jl")
    include("calculations/Conductance.jl")
end

# Run

key = ARGS[1]

if endswith(key, "_cond")
    truekey = replace(key, "_cond" => "")
    @info "Calculating conductance for system $(truekey)"
    res = calc_conductance(truekey)
    save(res.path, "res", res)
    @info "Conductance saved to $(res.path)"
    @info "Conductance calculation complete for system $(truekey)"
elseif endswith(key, "_spec")
    truekey = replace(key, "_spec" => "")
    @info "Calculating spectrum for system $(truekey)"
    res = calc_spectrum(truekey)
    save(res.path, "res", res)
    @info "Spectrum saved to $(res.path)"
    @info "Spectrum calculation complete for system $(truekey)"
else
    @error "Invalid key: $(key)"
end

# Epilog
rmprocs(workers()...)