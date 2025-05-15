@with_kw struct Lead_Params
    µ::Float64 = 0.0
    t::Float64 = 1
    τ::Float64 = 0.1
    dims::Int = 2
    nambu::Bool = true
end

function build_lead(params::Lead_Params)
    @unpack µ, t, τ, dims, nambu = params
    base = σ0
    if nambu 
        base = τz
    end
    if dims == 4
        base = σ0τz 
    end
    lat = LP.linear(; )
    model = onsite(µ * base) + hopping(t * base)
    h = lat |> hamiltonian(model; orbitals = dims)
    return h, τ
end