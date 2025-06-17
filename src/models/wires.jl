systems = Dict("base_wire" => System(;
        chain_params = Wire_Params(;
            α = 10.0,
            μ = 0.0,
            L = 500
            ),
        lead_params = Lead_Params(;
            t = 1.0,
            τ = 1.0
            ),
        params = Params(;
            Vzrng = range(0, .3, 501),
            ωrng = range(-.25, 0.4, 501) .+ 1e-4im,
            nev = 50,
            )
    )
)

systems["base_wire_nh"] = System(systems["base_wire"];
    NH_params = build_NH_params(
        systems["base_wire"].chain_params,
        [0.05, 0.05]
    )
)