@with_kw struct Wire_Params #meV, nm
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

build(params::Wire_Params) = build_wire(params)

function build_wire(params::Wire_Params)
    @unpack a0, t, N, α, µ, Vz, θ = params

    lat = LP.linear(; a0) |> supercell(N) |> supercell() 

    p2 = @onsite((r; µ = µ) -> (2 * t - μ) * σ0) + hopping((r, dr) -> -t * σ0)

    zeeman = @onsite((; Vz = Vz, θ = θ) -> Vz * cos(θ) * σx + Vz * sin(θ) * σy)

    rashba = @hopping((r, dr; α = α) -> α * im * σy * dr[1] / (2 * a0^2);)

    h = lat |> hamiltonian(p2 + zeeman + rashba; orbitals = 2)

    return h
end

@with_kw struct Wire_NH_Params
    γdict::Dict{Int, Vector{Float64}} = Dict(1 => [0.0, 0.0])
end

function build_NH_params(hparams::Wire_Params, γs)
    @unpack N = hparams
    if γs isa Dict{Int, Vector{Float64}}
        γdict = Dict(n => [0.0, 0.0] for n in 1:N)
        for (n, γ) in γs
            if n > N
                throw(ArgumentError("Key $n in γs is greater than the number of sites $N."))
            end
            γdict[n] = γ
        end
    elseif γs isa Vector{Float64} && length(γs) == 2
        γdict = Dict(n => copy(γs) for n in 1:N)
    else
        throw(ArgumentError("γs must be a dict or a vector of length 2."))
    end
    return Wire_NH_Params(; γdict)
end

add_NH_lead!(params::Wire_NH_Params) = h -> add_NH_lead!(h, params)
function add_NH_lead!(h::Quantica.AbstractHamiltonian, params::Wire_NH_Params)
    @unpack γdict = params
    
    sγ = γdict |> values |> sum
    if sγ == zeros(length(sγ))
        return h
    end
    sites = getindex.(h.hparent.lattice.unitcell.sites, 1)
    a0 = sites[2] - sites[1]

    if length(γdict) != length(sites)
        throw(ArgumentError("Length of γs must match the number of sites in the chain."))
    end

    sindex(r) = round(Int, r[1] / a0 + 1)
    nh_term = @onsite!((o, r;) -> o - 1im * (γdict[sindex(r)][1] * σ0 + γdict[sindex(r)][2] * σy); region = r -> (r[1] >= sites[1]) && (r[1] <= sites[end]))
    return h |> nh_term
end