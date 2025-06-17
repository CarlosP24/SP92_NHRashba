function calc_LDOS(name::String)
    system = systems[name]
    @unpack params, chain_params, NH_params = system
    @unpack ωrng, outdir, Vzrng = params

    h = build(chain_params) |> add_NH_lead!(NH_params)
    g = h |> greenfunction()

    ρup = ldos(g[], kernel = (σ0 + σy)/2)
    ρdown = ldos(g[], kernel = (σ0 - σy)/2)

    ρupf = (ω, Vz) -> ρup(ω; Vz)
    ρdownf = (ω, Vz) -> ρdown(ω; Vz)


    LDOSup = pfunction(ρupf, [ωrng, Vzrng])
    LDOSdown = pfunction(ρdownf, [ωrng, Vzrng])

    LDOS = (LDOSup, LDOSdown)
    path = "$(outdir)/LDOS/$(name).jld2"
    return Results(; system = system, LDOS = LDOS, path = path)
end