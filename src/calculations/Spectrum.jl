function calc_spectrum(name::String)
    system = systems[name]
    @unpack params, chain_params, NH_params= system
    if chain_params.N == 0
        throw(ArgumentError("chain_params.N cannot be 0"))
    end
    @unpack Vzrng, outdir, nev  = params

    h = build(chain_params) |> add_NH_lead!(NH_params)
    hf = Vz -> h(Vz = Vz)[()]

    Es = pspectrum(hf, [Vzrng...]; nev)

    path = "$(outdir)/Spectrum/$(name).jld2"
    return Results(; system = system, Es = Es, path = path)
end