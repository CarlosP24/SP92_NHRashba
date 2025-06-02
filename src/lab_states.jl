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


###
system = systems["Filter_real_nh_0.05"]

 @unpack params, chain_params, NH_params= system
if chain_params.N == 0
    throw(ArgumentError("chain_params.N cannot be 0"))
end
@unpack x, outdir, nev  = params

h = build(chain_params) |> add_NH_lead!(NH_params)

if x == :µ
    hf = µ -> h(µ = µ)[()]
    xrng = params.μrng
elseif x == :Vz
    hf = Vz -> h(Vz = Vz)[()]
    xrng = params.Vzrng
elseif x == :θ
    hf = θ -> h(θ = θ)[()]
    xrng = params.θrng
elseif x == (:Vz, :θ)
    hf = (Vz, θ) -> h(Vz = Vz, θ = θ)[()]
    xrng = (params.Vzrng, params.θrng)
else
    throw(ArgumentError("x must be :µ or :Vz"))
end

θ = 0
Vz = 0.15
λ, ψ = eigs(hf(Vz, θ); nev = 2000, which = :LM,sigma = 0.0, check = 1, tol = 1e-10,)
xrng = collect(1:chain_params.N).*chain_params.a0
ψodd = ψ[1:2:end, :]
ψeven = ψ[2:2:end, :]
λodd = λ[1:2:end]
λeven = λ[2:2:end]

##
ψp = 0.5*sqrt.(abs.(ψodd).^2 + abs.(ψeven).^2)

##
# Sort λ and ψp_mat columns by real(λ)
idx = sortperm(real.(λ))
λs = λ[idx]
ψs = ψp[:, idx] ./ maximum(ψp)

sample = 10
λp = λs[1:sample:end]
ψpp = ψs[:, 1:sample:end]

##
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"$x$ (nm)", ylabel = L"$|\psi|$")
colors = cgrad(:redsblues, length(λp), categorical = true)

lbig = round(Int ,length(λp) * 0.01)
αs = vcat(ones(lbig), 0.05 * ones(length(λp) - 2*lbig), ones(lbig))

sorted_ψpps = hcat(ψpp[:, 1:lbig], ψpp[:, (end-lbig+1):end], ψpp[:, (lbig+1):(end-lbig)])
sorted_λs = real.(vcat(λp[1:lbig], λp[(end-lbig+1):end], λp[(lbig+1):(end-lbig)]))
mλ = maximum(sorted_λs)

for (i, λ) in enumerate(sorted_λs)
    lines!(ax, xrng, sorted_ψpps[:, i]; color = (colors[ceil(Int, length(λp) *    λ/mλ|>abs)], αs[ceil(Int,length(λp) * λ/mλ |>abs)]), linewidth = 2)
end 
fig

##
t = ħ2ome / (2 * m0 * a0^2)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"$\text{Re}(\omega - 2t)/(2t)$", ylabel = L"$\text{Im}(\omega)$ (meV)")
scatter!(ax, (real.(λ) .- 2*t)/(2*t), imag.(λ); color = :black)
fig



##
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"$x$ (nm)", ylabel = L"$\text{Re}(\omega)$ (meV)")
heatmap!(ax, xrngp, real.(λp), ψpp; colormap = :thermal, colorrange = (0, 1))
ϵ1 = 0
iϵ1 = findmin(abs.(real.(λp) .- ϵ1))[2]

ϵ3 = 154
iϵ3 = findmin(abs.(real.(λp) .- ϵ3))[2]

ϵ2 = 266
iϵ2 = findmin(abs.(real.(λp) .- ϵ2))[2]

hlines!(ax, ϵ1; color = :red, linestyle = :dash)
hlines!(ax, ϵ2; color = :green, linestyle = :dash)
hlines!(ax, ϵ3; color = :blue, linestyle = :dash)
hidexdecorations!(ax, ticks = false)

ax = Axis(fig[2, 1]; xlabel = L"$x$ (nm)", ylabel = L"$|\psi|$")
lines!(ax, xrngp, ψpp[:, iϵ1]; color = :red)
lines!(ax, xrngp, ψpp[:, iϵ2]; color = :green)
lines!(ax, xrngp, ψpp[:, iϵ3]; color = (:blue, 0.1))
xlims!(ax, (first(xrngp), last(xrngp)))
#ylims!(ax, (minimum(real.(λ)), 10))
#hlines!(ax, 0.05; color = :white, linestyle = :dash)
Colorbar(fig[1, 2], colormap = :thermal, limits = (0, 1), ticks = [0, 1], label = L"$|\psi|$", labelpadding = -15)

# ax = Axis(fig[2, 1]; xlabel = L"$\text{Re}(\omega)$ (meV)", ylabel = L"$\text{Im}(\omega)$ (meV)")
# scatter!(ax, real.(λ), imag.(λ); color = :black)
rowgap!(fig.layout, 1, 5)
colgap!(fig.layout, 1, 5)

fig