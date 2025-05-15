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
    nambu::Bool = false
end

build(params::Filter_Params) = build_filter(params)

function build_filter(params::Filter_Params)
    @unpack a0, t, N, α, µ, Vz, θ = params

    lat = LP.linear(; a0) |> supercell(N) |> supercell() 

    p2 = @onsite((r; µ = µ) -> (2 * t - μ) * σ0) + hopping((r, dr) -> -t * σ0)

    zeeman = @onsite((; Vz = Vz, θ = θ) -> Vz * cos(θ) * σx + Vz * sin(θ) * σy)

    rashba = @hopping((r, dr; α = α) -> α * im * σy * dr[1] / (2 * a0^2);)

    h = lat |> hamiltonian(p2 + zeeman + rashba; orbitals = 2)

    return h
end

@with_kw struct Filter_NH_Params
    γ0s::Vector{Float64} = [0.0]
    γys::Vector{Float64} = [0.0]
end

function build_NH_params(hparams::Filter_Params, γs::Dict{Int, Vector{Float64}})
    @unpack N = hparams
    γ0s = zeros(Float64, N)
    γys = zeros(Float64, N)
    for (i, γs) in γs
        γ0s[i] = γs[1]
        γys[i] = γs[2]
    end
    return Filter_NH_Params(γ0s, γys)
end

add_NH_lead!(params::Filter_NH_Params) = h -> add_NH_lead!(h, params)
function add_NH_lead!(h::Quantica.AbstractHamiltonian, params::Filter_NH_Params)
    @unpack γ0s, γys = params
    if γ0s == [0.0] && γys == [0.0]
        return h
    end
    sites = getindex.(h.hparent.lattice.unitcell.sites, 1)

    if length(γ0s) != length(sites)
        throw(ArgumentError("Length of γ0s must match the number of sites in the chain."))
    end
    if length(γys) != length(sites)
        throw(ArgumentError("Length of γys must match the number of sites in the chain."))
    end

    for (γ0, γy, site) in zip(γ0s, γys, sites)
        if γ0 == 0 && γy == 0
            continue
        end
        nh_term = @onsite!((o, r;) -> o - 1im * (γ0 * σ0 + γy * σy); region = r -> r[1] == site)
        h = h |> nh_term
    end
    return h
end