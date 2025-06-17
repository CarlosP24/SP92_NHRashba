function calc_conductance(name::String)
    system = systems[name]
    @unpack params, chain_params, NH_params, lead_params = system
        if chain_params.N == 0
        throw(ArgumentError("chain_params.N cannot be 0"))
    end
    @unpack a0 = chain_params
    @unpack ωrng, outdir, Vzrng = params

    h = build(chain_params) |> add_NH_lead!(NH_params)
    h_lead, τt = build_lead(lead_params)

    glead = h_lead |> greenfunction(GS.Schur(boundary = 0))
    g = h |> attach(glead, @hopping((; τt = τt)-> τt * σ0), region = r -> r[1] == 0) |> attach(glead, @hopping((; τt = τt)-> τt * σ0), region = r -> r[1] == (chain_params.N-1) * a0) |> greenfunction()

    gpts = Iterators.product(1:2, 1:2)

    Gs = map(gpts) do (i, j)
        G = conductance(g[i, j];)
        Gf = (ω, Vz) -> G(ω; Vz)
        return pfunction(Gf, [ωrng, Vzrng])
    end

    Gs = reshape(Gs, 2, 2)

    path = "$(outdir)/Conductance/$(name).jld2"
    return Results(; system = system, Gs = Gs, path = path)
end