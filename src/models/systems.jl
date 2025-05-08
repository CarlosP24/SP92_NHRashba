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
    chain_params::Kitaev_Params = Kitaev_Params(; )
    NH_params::NH_Lead_Params = NH_Lead_Params(; )
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