@with_kw struct Lead_Params
    µ::Float64 = 0.0
    t::Float64 = 1
    τ::Float64 = 0.1
end

function build_lead(params::Lead_Params)
    @unpack µ, t, τ = params
    lat = LP.linear(; )
    model = onsite(µ * τz) + hopping(t * τz)
    h = lat |> hamiltonian(model; orbitals = 2)
    return h, τ
end