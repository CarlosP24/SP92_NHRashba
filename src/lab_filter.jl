using Quantica, Parameters
using Quantica: σ
@with_kw struct Filter_Params #meV, nm
    α::Float64 = 0.0
    μ::Float64 = 0.0
    m0::Float64 = 0.023
    a0::Float64 = 5
    ħ2ome::Float64 = 76.1996
    L::Float64 = 100
    N::Int = round(Int, L / a0)
    t = ħ2ome / (2 * m0 * a0^2)
    Vz::Float64 = 0.0
    θ::Float64 = 0.0
end


σz = σ(3)
σx = σ(1)
σy = σ(2)
σ0 = σ(0)

function build_filter(params::Filter_Params)
    @unpack a0, t, N, α, µ, Vz, θ = params

    lat = LP.linear(; a0) |> supercell(N) |> supercell() 

    p2 = @onsite((r; µ = µ) -> (2 * t - μ) * σ0) + hopping((r, dr) -> -t * σ0)

    zeeman = @onsite((; Vz = Vz, θ = θ) -> Vz * cos(θ) * σx + Vz * sin(θ) * σy)

    rashba = @hopping((r, dr; α = α) -> α * im * σy * dr[1] / (2 * a0^2);)

    h = lat |> hamiltonian(p2 + zeeman + rashba; orbitals = 2)

    return h
end

function build_lead_normal(params::Lead_Params)
    @unpack µ, t, τ  = params
    lat = LP.linear(; )
    model = onsite(µ * σ0) + hopping(t * σ0)
    h = lat |> hamiltonian(model; orbitals = 2)
    return h, τ
end

h_lead, τ = build_lead_normal(Lead_Params(; ))
fparams = Filter_Params(; )
h = build_filter(fparams)
glead = h_lead |> greenfunction(GS.Schur(boundary = 0))

g = h |> attach(glead, @hopping((; τ = τ) -> τ * σ0), region = r -> r[1] == 0) |> attach(glead, @hopping((; τ = τ) -> τ * σ0), region = r -> r[1] == (fparams.N-1) * fparams.a0) |> greenfunction()