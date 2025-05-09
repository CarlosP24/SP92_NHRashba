@with_kw struct Params
    µrng = subdiv(-4, 4, 401)
    ωrng = subdiv(-5, 5, 401) .+ 1e-2im
    outdir = "data"
end

@with_kw struct Results
    system = nothing
    Gs = nothing
    Es = nothing
    path = nothing
end

@with_kw struct System
    params::Params = Params()
    chain_params = Kitaev_Params(; )
    NH_params = NH_Lead_Params(; )
    lead_params::Lead_Params = Lead_Params(; )
end

systems = Dict(
    "hermitian_chain" => System(
        chain_params = Kitaev_Params(; N = 50, Δ = 2),
        NH_params = NH_Lead_Params(; Γodd = 0.0, Γeven = 0.0),
        lead_params = Lead_Params(; ),
    ),
)

systems["nh_even"] = System(systems["hermitian_chain"]; 
    NH_params = NH_Lead_Params(; Γodd = 0.0, Γeven = 0.1),
)

systems["nh_odd"] = System(systems["hermitian_chain"]; 
    NH_params = NH_Lead_Params(; Γodd = 1.0, Γeven = 0.0),
)

systems["nh_odd_left"] = System(systems["hermitian_chain"]; 
    NH_params = NH_Lead_Params(; Γodd_left = 1.1, Γodd_right = 0.9, Γeven = 0.0),
)

systems["nh_odd_right"] = System(systems["hermitian_chain"]; 
    NH_params = NH_Lead_Params(; Γodd_left = 0.9, Γodd_right = 1.1, Γeven = 0.0),
)

systems["nh_odd_superleft"] = System(systems["hermitian_chain"]; 
    NH_params = NH_Lead_Params(; Γodd_left = 2, Γodd_right = 1, Γeven = 0.0),
)

params_wire = Params(;
    ωrng = subdiv(-.25, .25, 401) .+ 1e-3im,
    μrng = subdiv(0, 5, 401),  
)

systems["Wire_base"] = System(
    chain_params = Wire_Params(; L = 2000, Δ = 0.23),
    NH_params = wire_nh_lead_params(; γ0 = 0.0, γy = 0.0),
    params = params_wire,
)