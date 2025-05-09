
@with_kw struct Wire_Params # units meV, nm, T
    ħ2ome::Float64 = 76.1996
    #μB::Float64 = 5.78e-2               
    m0::Float64 = 0.023
    #g::Float64 = 10
    α::Float64 = 5
    μ::Float64 = 0.0
    a0::Float64 = 5
    t = ħ2ome / (2 * m0 * a0^2)
    Δ::Float64 = 0.23
    L::Float64 = 0
    N = round(Int, L / a0)
    #B = 1
    Vz = 0.0
end

build(params::Wire_Params) = build_wire(params)

function build_wire(p::Wire_Params)
    @unpack ħ2ome, m0, α, µ, a0, t, Δ, N, Vz = p

    lat = LP.linear(; a0)
    
    if N != 0 
        lat = lat |> supercell(N) |> supercell()
    end

    p2 = @onsite((r; µ = µ) -> (2 * t - µ) * σ0τz ) + hopping((r, dr; t = t) -> -t * σ0τz; range = a0)

    zeeman = @onsite((; Vz = Vz) -> σzτ0 * Vz)

    rashba = @hopping((r, dr; α = α) -> σyτz * α * im * dr[1] / (2 * a0^2); range = a0)

    pairing = @onsite((; Δ = Δ) -> Δ * σ0τx)

    h = lat |> hamiltonian(p2 + zeeman + rashba + pairing; orbitals = 4)

    return h
end

@with_kw struct wire_nh_lead_params
    γ0::Float64 = 0.0
    γx::Float64 = 0.0
    γy::Float64 = 0.0
end

add_NH_lead!(params::wire_nh_lead_params) = h -> add_NH_lead!(h, params)
function add_NH_lead!(h::Quantica.AbstractHamiltonian, params::wire_nh_lead_params)
    @unpack γ0, γx, γy = params

    nh_term = @onsite!((o, r; γ0 = γ0, γy = γy) -> o + 1im * (γ0 * σ0τ0 + γy * σyτ0 + γx * σxτ0); )

    return h |> nh_term
end