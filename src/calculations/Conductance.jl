function calc_conductance(name::String)
    system = systems[name]
    @unpack params, chain_params, NH_params, lead_params = system
    @unpack ωrng, μrng, outdir = params

    h = build_Kitaev(chain_params) |> add_NH_lead!(NH_params)
    h_lead, τ = build_lead(lead_params)

    glead = h_lead |> greenfunction(GS.Schur(boundary = 0))
    g = h |> attach(glead, @hopping((; τ = τ)-> τ * τz), region = r -> r[1] == 0) |> attach(glead, @hopping((; τ = τ)-> τ * τz), region = r -> r[1] == (chain_params.N-1)) |> greenfunction()

    gpts = Iterators.product(1:2, 1:2)

    Gs = map(gpts) do (i, j)
        G = conductance(g[i, j]; nambu = true)
        return pfunction((ω, μ) -> G(ω; μ), [ωrng, μrng])
    end

    Gs = reshape(Gs, 2, 2)

    path = "$(outdir)/Conductance/$(name).jld2"
    return Results(; system = system, Gs = Gs, path = path)
end