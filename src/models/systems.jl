@with_kw struct Params
    Vzrng = range(0, 4, 201)
    ωrng = range(-5, 5, 201) .+ 1e-2im
    outdir = "data"
    nev = 20
end
@with_kw struct System
    params::Params = Params()
    chain_params = Wire_Params(; )
    NH_params = Wire_NH_Params(; )
    lead_params::Lead_Params = Lead_Params(; )
end
@with_kw struct Results
    system = nothing
    Gs = nothing
    Es = nothing
    LDOS = nothing
    path = nothing
end

include("wires.jl")