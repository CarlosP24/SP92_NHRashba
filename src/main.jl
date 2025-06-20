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

    include("builders/Wire.jl")
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
elseif key in keys(systems)
    @info "Calculating spectrum and conductance for system $(key)"
    res = calc_spectrum(key)
    save(res.path, "res", res)
    @info "Spectrum saved to $(res.path)"
    @info "Conductance calculation for system $(key)"
    res = calc_conductance(key)
    save(res.path, "res", res)
    @info "Conductance saved to $(res.path)"
    @info "Spectrum and conductance calculation complete for system $(key)"
    @info "LDOS calculation for system $(key)"
    res = calc_LDOS(key)
    save(res.path, "res", res)
    @info "LDOS saved to $(res.path)"
else
    @error "Invalid key: $(key)"
end

# Epilog
rmprocs(workers()...)