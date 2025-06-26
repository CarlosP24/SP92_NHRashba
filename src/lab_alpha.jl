using Distributed
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

cd("src")

include("builders/Wire.jl")
include("builders/leads.jl")

include("models/systems.jl")
include("parallelizers/pgeneric.jl")
include("parallelizers/pspectrum.jl")

include("calculations/Spectrum.jl")
include("calculations/Conductance.jl")
#include("calculations/LDOS.jl")

##
name = "base_wire_nh"
system = systems[name]
@unpack params, chain_params, NH_params = system

h = build(chain_params) |> add_NH_lead!(NH_params)
hf = (Vz, α) -> h(Vz = Vz, α = α)[()]

Vzrng = range(0, 0.3, 201)
αrng = range(0.0, 20.0, 201)
nev = 6
Es = pspectrum(hf, [Vzrng...],[αrng...]; nev)

##
αs = [0, 10, 20]
colors = [:red, :blue, :green]
fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$V_Z$ (meV)", ylabel = L"$E$ (meV)")
for n in 1:nev
    for (α, color) in zip(αs, colors)
        iα = findmin(abs.(αrng .- α))[2]
        E = Es[:, iα]
        scatter!(ax, Vzrng, real.(getindex.(E, n)), color = color)
    end
end
fig
##
fig = Figure()
for (n1, n2) in zip(1:2:nev, 2:2:nev)
    ΔImEs = [imag(e[n2]) - imag(e[n1]) for e in Es]
    ax = Axis(fig[Int(ceil(n1/2)), 1], xlabel = L"$V_Z$ (meV)", ylabel = L"$\alpha$ (meVnm)")
    heatmap!(ax, Vzrng, αrng, abs.(ΔImEs), colormap = :amp)
    contour!(ax, Vzrng, αrng, abs.(ΔImEs), levels = [1e-3])
    vlines!(ax, 0.05; color = :white, linestyle = :dash)
    n2 != nev && hidexdecorations!(ax, ticks = false, grid = false)
    n1 != 1 && rowgap!(fig.layout, Int(ceil(n1/2))-1, 5)
    ylims!(ax, 0, last(αrng))
    xlims!(ax, 0, last(Vzrng))
end
fig
