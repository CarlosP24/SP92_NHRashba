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
using CairoMakie


include("builders/Kitaev.jl")
include("builders/Wire.jl")
include("builders/Filter.jl")
include("builders/leads.jl")

include("models/systems.jl")
include("models/Filters.jl")
###
function sumψ(ψ)
    ψodd = ψ[1:2:end]
    ψeven = ψ[2:2:end]
    return 0.5*sqrt.(abs.(ψodd).^2 + abs.(ψeven).^2) |> vec
end
    ###
name = "F_N=100"
system = systems[name]
@unpack params, chain_params, NH_params= system

h = build(chain_params) |> add_NH_lead!(NH_params)
hf = (Vz) -> h(Vz = Vz)[()]

###
Vzrng = range(0, 0.3, 201)
nev = 4
EsR = Dict{Float64, Vector{ComplexF64}}()
ψsR = Dict{Float64, Matrix{ComplexF64}}()
EsL = Dict{Float64, Vector{ComplexF64}}()
ψsL = Dict{Float64, Matrix{ComplexF64}}()
for Vz in Vzrng
    λR, ψR = eigs(hf(Vz); nev, ncv = nev, which = :LM, sigma = 0.0, check = 1, tol = 1e-10)
    λL, ψL = eigs(hf(Vz)'; nev, ncv = nev, which = :LM, sigma = 0.0, check = 1, tol = 1e-10)
    idxL = sortperm(real(λL))
    λL = λL[idxL]
    ψL = ψL[:, idxL]
    idxR = sortperm(real(λR))
    λR = λR[idxR]
    ψR = ψR[:, idxR]
    EsR[Vz] = λR
    ψsR[Vz] = ψR
    EsL[Vz] = conj.(λL)
    ψsL[Vz] = conj.(ψL)
end

###
for (i, Vz) in enumerate(Vzrng)
    Vz_closest = Vzrng[argmin(abs.(Vzrng .- Vz))]
    fig = Figure(size = (1000, 500), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = L"$V_Z$ (meV)", ylabel = L"$E$ (meV)")
    for Vz in Vzrng, n in 1:4
        scatter!(ax, Vz, real(EsR[Vz][n]), color = :green)
        scatter!(ax, Vz, imag(EsR[Vz][n]), color = :darkred, markersize = 2)
    end
    Label(fig[1, 1, Top()], L"$V_Z = %$(round(Vz_closest, digits = 2)) \text{ meV}$", padding = (-250, 0, -50, 0), fontsize = 20)

    fig_ψs = fig[1, 2] = GridLayout()
    xlims!(ax, first(Vzrng), last(Vzrng))
    vlines!(ax, Vz_closest; color = :black, linestyle = :dash)
    xrng = collect(1:chain_params.N) .* chain_params.a0

    for n in 1:nev
        ax = Axis(fig_ψs[nev-n + 1, 1], xlabel = L"$x$ (nm)", ylabel = L"$|\psi|$", yticks = [0])
        ψR = ψsR[Vz_closest][:, n]
        ψL = ψsL[Vz_closest][:, n]
        lines!(ax, xrng, sumψ(ψR), color = :blue, label = "Right" )
        lines!(ax, xrng, sumψ(ψL), color = :red, label = "Left")
        xlims!(ax, first(xrng), last(xrng))
        ylims!(ax, 0, 1.1*maximum([sumψ(ψR); sumψ(ψL)]))
        n == 1 && axislegend(ax; position = :cb, orientation = :horizontal, framevisible = false)
        n != 1 && hidexdecorations!(ax, ticks = false, grid = false)
        Label(fig_ψs[nev-n + 1, 1, Top()], L"$n = %$(n)$", padding = (-480, 0, -20, 0), fontsize = 20)
    end
    save("plots/figures/F_N=100_Vzfilm/$(i)_Vz=$(Vz).pdf", fig)
end