function calc_spectrum(name::String)
    system = systems[name]
    @unpack params, chain_params, NH_params= system
    if chain_params.N == 0
        throw(ArgumentError("chain_params.N cannot be 0"))
    end
    @unpack x, outdir, nev  = params

    h = build(chain_params) |> add_NH_lead!(NH_params)

    if x == :µ
        hf = µ -> h(µ = µ)[()]
        xrng = params.μrng
    elseif x == :Vz
        hf = Vz -> h(Vz = Vz)[()]
        xrng = params.Vzrng
    elseif x == :θ
        hf = θ -> h(θ = θ)[()]
        xrng = params.θrng
    elseif x == (:Vz, :θ)
        hf = (Vz, θ) -> h(Vz = Vz, θ = θ)[()]
        xrng = (params.Vzrng, params.θrng)
    else
        throw(ArgumentError("x must be :µ or :Vz"))
    end
    
    if length(xrng) == 1
        Es = pspectrum(hf, xrng; nev)
    else
        Es = pspectrum(hf, xrng...; nev)
    end

    path = "$(outdir)/Spectrum/$(name).jld2"
    return Results(; system = system, Es = Es, path = path)
end