@with_kw struct Lead_Params
    µ::Float64 = 0.0
    t::Float64 = 1
    τ::Float64 = 0.1
    dims::Int = 2
end

function build_lead(params::Lead_Params)
    @unpack µ, t, τ, dims = params
    if dims == 4
        τz = σ0τz 
    end
    lat = LP.linear(; )
    model = onsite(µ * τz) + hopping(t * τz)
    h = lat |> hamiltonian(model; orbitals = dims)
    return h, τ
end