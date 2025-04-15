using Quantica, Parameters, CairoMakie
using LinearAlgebra, Arpack

const τ0 = @SMatrix[1 0; 0 1]
const τx = @SMatrix[0 1; 1 0]
const τy = @SMatrix[0 -1im; 1im 0]
const τz = @SMatrix[1 0; 0 -1]

@with_kw struct Kitaev_Params
    μ::Float64 = 0.0
    t::Float64 = 1.0
    Δ::Float64 = 1.0 
    N::Int = 25
    Γ::Float64 = 0.0
end

function build_Kitaev(params::Kitaev_Params)
    @unpack µ, t, Δ, N, Γ = params

    lat = LP.linear(; ) |> supercell(N) |> supercell()

    kinetic = @onsite((; µ = 0) -> -µ * τz) + @hopping((; t = t) -> -t * τz)
    pairing1 = @hopping((; Δ = Δ) -> -Δ * 1im * τy;
    region = (r, dr) -> dr[1] > 0 )
    pairing2 = @hopping((; Δ = Δ) -> Δ * 1im * τy;
    region = (r, dr) -> dr[1] < 0 )
    nh_h = @hopping((r, dr; Γ = Γ) -> Γ * τ0; region = (r, dr) -> dr[1] > 0) + 
    @hopping((r, dr; Γ = Γ) -> - Γ * τ0; region = (r, dr) -> dr[1] < 0)

    h = lat |> hamiltonian(kinetic + pairing1 + pairing2 + nh_h; orbitals = 2)
    return h 
end

@with_kw struct Lead_Params
    µ::Float64 = 0.0
    t::Float64 = 1.0
    τ::Float64 = 0.1
end

function build_lead(params::Lead_Params)
    @unpack µ, t, τ = params
    lat = LP.linear(; )
    model = onsite(-µ * τz) + hopping(-t * τz)
    h = lat |> hamiltonian(model; orbitals = 2)
    return h, τ
end

function get_spectrum(hfunc, µrng; nev = 20)
    Es = map(µrng) do µ
        try
            λ, ψ = eigs(hfunc(µ); nev, ncv = nev , which = :LM, sigma = 0.0, check = 1, tol = 1e-10)
        catch
            return fill(NaN, nev)
        end
        if length(λ) < nev
            λ = vcat(λ, fill(NaN, nev - length(λ)))
        end
        return λ
    end
    Es = reshape(Es, size(µrng)...)
    Es = hcat(Es...)
    sort!(Es, dims = 1, by = x -> abs(real(x)))
    return Es
end

function get_states(h; nev = 20)
    λ, ψ = eigs(h[()]; nev, ncv = nev , which = :LM, sigma = 0.0, check = 1)
    ψup = ψ[1:2:end, :]
    ψdown = ψ[2:2:end, :]
    return ψup, ψdown
end

function nor_sum(u, v)
    s = abs.(u).^2 + abs.(v).^2
    maxs = maximum(s)
    return vec(s)
end
##
kchain = Kitaev_Params(; N = 25, Δ = 2)
h = build_Kitaev(kchain)

µrng = subdiv(0, 4, 201)

Es = get_spectrum(µ -> h(µ = µ)[()], µrng; nev = 20)
Es_nh = get_spectrum(µ -> h(µ = µ, Γ = 1)[()], µrng; nev = 20)
ψup, ψdown = get_states(h(µ = 1.0); )
ψup_nh, ψdown_nh = get_states(h(µ = 1.0, Γ = 1); )


## Conductance
hlead, τ = build_lead(Lead_Params(;))
glead = hlead |> greenfunction(GS.Schur(boundary = 0))
g = h |> attach(glead, @hopping((; τ = τ)-> τ * τz), region = r -> r[1] == 0) |> attach(glead, @hopping((; τ = τ)-> τ * τz), region = r -> r[1] == (kchain.N-1)) |> greenfunction()

GLL = conductance(g[1, 1])
GRR = conductance(g[2, 2])
GLR = conductance(g[1, 2])
GRL = conductance(g[2, 1])

function get_conductance(Gfunc, ωrng, µrng)
    pts = Iterators.product(µrng, ωrng)
    Gs = map(pts) do (μ, ω)
        return Gfunc(ω, μ)
    end
    return Gs
end
ωrng = subdiv(-1, 1, 201) .+ 1e-3im
GLLs = get_conductance((ω, µ) -> GLL(ω; µ), ωrng, µrng)
GRRs = get_conductance((ω, µ) -> GRR(ω; µ), ωrng, µrng)
GRLs = get_conductance((ω, µ) -> GRL(ω; µ), ωrng, µrng)
GLRs = get_conductance((ω, µ) -> GLR(ω; µ), ωrng, µrng)

GLLs_nh = get_conductance((ω, µ) -> GLL(ω; µ, Γ = 1), ωrng, µrng)
GRRs_nh = get_conductance((ω, µ) -> GRR(ω; µ, Γ = 1), ωrng, µrng)
GRLs_nh = get_conductance((ω, µ) -> GRL(ω; µ, Γ = 1), ωrng, µrng)
GLRs_nh = get_conductance((ω, µ) -> GLR(ω; µ, Γ = 1), ωrng, µrng)


## Plotting
function plot_spectrum(pos, µrng, Es)
    ax = Axis(pos; xlabel = L"\mu / t", ylabel = L"E / t")
    color = :blue
    color_imag = :purple
    for (i, E) in enumerate(eachrow(Es))
        scatter!(ax, µrng, real.(E); color, markersize = 5)
        scatter!(ax, µrng, imag.(E); color = color_imag, markersize = 2)
        if i == 2
            color = :black
            color_imag = :red
        end
    end
    ylims!(ax, (-2, 2))
    xlims!(ax, (0, 4))
    return ax
end

function plot_state(pos, µrng, ψup, ψdown)
    ax = Axis(pos; xlabel = L"n", ylabel = L"|u|^2 + |v|^2")
    lines!(ax, nor_sum(ψup, ψdown); color = :darkgreen, linewidth = 3)
    return ax
end
fig = Figure()
ax = plot_spectrum(fig[1, 1], µrng, Es)
vlines!(ax, 1; color = :darkgreen, linestyle = :dash)
hidexdecorations!(ax, ticks = false, grid = false)
ax = plot_state(fig[1, 2], µrng, ψup, ψdown)
xlims!(ax, (1, 25))
hidexdecorations!(ax, ticks = false, grid = false)

ax = plot_spectrum(fig[2, 1], µrng, Es_nh)
vlines!(ax, 1; color = :darkgreen, linestyle = :dash)

ax = plot_state(fig[2, 2], µrng, ψup_nh[:, 1], ψdown_nh[:, 1])
xlims!(ax, (1, 25))
fig
##
function plot_conductance(pos, Gs, ωrng, μrng; colorrange = (-.1, .1))
    ax = Axis(pos; xlabel = L"\mu / t", ylabel = L"\omega / t")
   hmap = heatmap!(ax, μrng, real.(ωrng), Gs; colormap = :balance, colorrange)
    return ax, hmap
end

fig = Figure()
ax, hmap = plot_conductance(fig[1, 1], GLLs, ωrng, µrng)
hidexdecorations!(ax, ticks = false, grid = false)
ax, hmap = plot_conductance(fig[1, 2], GLRs, ωrng, µrng)
hidexdecorations!(ax, ticks = false, grid = false)
hideydecorations!(ax, ticks = false, grid = false)
ax, hmap = plot_conductance(fig[2, 1], GRLs, ωrng, µrng)
ax, hmap = plot_conductance(fig[2, 2], GRRs, ωrng, µrng)
hideydecorations!(ax, ticks = false, grid = false)
fig

##
fig = Figure()
ax, hmap = plot_conductance(fig[1, 1], GLLs_nh, ωrng, µrng)
hidexdecorations!(ax, ticks = false, grid = false)
ax, hmap = plot_conductance(fig[1, 2], GLRs_nh, ωrng, µrng)
hidexdecorations!(ax, ticks = false, grid = false)
hideydecorations!(ax, ticks = false, grid = false)
ax, hmap = plot_conductance(fig[2, 1], GRLs_nh, ωrng, µrng)
ax, hmap = plot_conductance(fig[2, 2], GRRs_nh, ωrng, µrng)
hideydecorations!(ax, ticks = false, grid = false)
fig
