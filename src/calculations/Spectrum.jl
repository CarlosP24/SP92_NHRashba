function calc_spectrum(name::String)
    system = systems[name]
    @unpack params, chain_params, NH_params = system
    if chain_params.N == 0
        throw(ArgumentError("chain_params.N cannot be 0"))
    end
    @unpack μrng, outdir = params

    h = build(chain_params) |> add_NH_lead!(NH_params)
    Es = pspectrum(µ -> h(µ = µ)[()], μrng; nev = 20)

    path = "$(outdir)/Spectrum/$(name).jld2"
    return Results(; system = system, Es = Es, path = path)
end