# Filter study loop
systems["F_Base"] = System(;
    chain_params = Filter_Params(; 
        α = 10.0,
        μ = 0.0,
    ),
    NH_params = Filter_NH_Params(;),
    params = Params(;
        ωrng = subdiv(-.25, .25, 801) .+ 1e-3im,
        Vzrng = subdiv(0, .2, 1001), 
        τrng = [1.0],
        x = (:Vz, :τ),
        nev = 50,
    ),
    lead_params = Lead_Params(; 
        nambu = false,
        t = 1)
)

Nrng = [10, 50, 100, 500, 1000]

for N in Nrng
    systems["F_N=$(N)"] = System(systems["F_Base"];
        chain_params = Filter_Params(systems["F_Base"].chain_params;
            N = N, 
            L = N * systems["F_Base"].chain_params.a0
        ),
        NH_params = Filter_NH_Params(systems["F_Base"].NH_params;
            γdict = Dict(i => [0.05, 0.05] for i in 1:N),
        ),
    )     
end