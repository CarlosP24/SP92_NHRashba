function calc_conductance(name::String)
    system = systems[name]
    @unpack params, chain_params, NH_params, lead_params = system
    if chain_params.N == 0
        throw(ArgumentError("chain_params.N cannot be 0"))
    end
    @unpack ωrng, outdir, x = params

    a0 = 1

    if chain_params isa Wire_Params
        lead_params = Lead_Params(lead_params; dims = 4)
        τz = σ0τz
        a0 = chain_params.a0
    end

    h = build(chain_params) |> add_NH_lead!(NH_params)
    h_lead, τ = build_lead(lead_params)

    glead = h_lead |> greenfunction(GS.Schur(boundary = 0))
    g = h |> attach(glead, @hopping((; τ = τ)-> τ * τz), region = r -> r[1] == 0) |> attach(glead, @hopping((; τ = τ)-> τ * τz), region = r -> r[1] == (chain_params.N-1) * a0) |> greenfunction()

    gpts = Iterators.product(1:2, 1:2)

    Gs = map(gpts) do (i, j)
        G = conductance(g[i, j]; nambu = true)
        if x == :µ
            Gf = (ω, μ) -> G(ω; μ)
            xrng = params.µrng
        elseif x == :Vz
            Gf = (ω, Vz) -> G(ω; Vz)
            xrng = params.Vzrng
        else
            throw(ArgumentError("x must be :µ or :Vz"))
        end
        return pfunction(Gf, [ωrng, xrng])
    end

    Gs = reshape(Gs, 2, 2)

    path = "$(outdir)/Conductance/$(name).jld2"
    return Results(; system = system, Gs = Gs, path = path)
end