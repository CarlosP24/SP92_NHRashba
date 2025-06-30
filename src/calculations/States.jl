    function calc_states(name::String, Vz::Float64, nev::Int)
        system = systems[name]
        @unpack params, chain_params, NH_params= system

        chain_params = Wire_Params(chain_params; Vz = Vz)
        systemV = System(system; chain_params = chain_params)
        @unpack outdir = params


        h = build(chain_params) |> add_NH_lead!(NH_params)
        hmat = h(;)[()]
        
        λR, ΨR = eigs(hmat; nev = nev, ncv = nev, which = :LM, sigma = 0.0, check = 1, tol = 1e-10,)
        idxR = sortperm(real(λR))
        λR = λR[idxR]
        ΨR = ΨR[:, idxR]

        λL, ΨL = eigs(hmat'; nev = nev, ncv = nev, which = :LM, sigma = 0.0, check = 1, tol = 1e-10,)
        idxL = sortperm(real(λL))
        λL = conj.(λL)
        λL = λL[idxL]
        ΨL = ΨL[:, idxL]

        ΨL = conj.(ΨL)

        path = "$(outdir)/States/$(name).jld2"
        return Results(; system = systemV, Ψs = [ΨR, ΨL], path = path )
    end