function fig_GVz(name::String; maxG = 1, τ = 1, enhance = 1)
    res = load("data/Conductance/$(name).jld2")["res"];
    eres = load("data/Spectrum/$(name).jld2")["res"];
    @unpack Es = eres;
    @unpack system, Gs, path = res;
    @unpack params = system;
    @unpack ωrng, Vzrng, τrng = params;

    γ = system.NH_params.γdict[1][2]

    τi = findmin(abs.(τrng .- τ))[2]
    fig = Figure()

    for i in axes(Gs, 2), j in axes(Gs, 1)
        G = Gs[i, j][:, :, τi]
        if i != j 
            G *= enhance
        end
        ax, hmap = plot_conductance(fig[i, j], G, ωrng, Vzrng; colorrange = (-maxG, maxG), labels = labels["Filter"])
        
        ax.xticks = ([first(Vzrng), γ, last(Vzrng)], [L"0", L"\gamma_y", L"%$(last(Vzrng))"])
        ax.xlabel = L"$V_Z$ (meV)"
        ax.yticks = [first(ωrng), 0, last(ωrng)] .|> real
        ax.ylabel = L"$\omega$ (meV)"
        vlines!(ax, system.NH_params.γdict[1][2]; color = :black, linestyle = :dash)
        ylims!(ax, (first(ωrng) |> real, last(ωrng) |> real))
        xlims!(ax, (first(Vzrng), last(Vzrng)))

        j == 2 && hideydecorations!(ax, ticks = false, grid = false)
        i == 1 && hidexdecorations!(ax, ticks = false, grid = false)

        i == j && Label(fig[i, j, Top()], L"$G_{%$(contact_dict[i]) %$(contact_dict[j])}$", fontsize = 15, padding = (180, 0, -300, 0))
        i != j && Label(fig[i, j, Top()], L"$G_{%$(contact_dict[i]) %$(contact_dict[j])} \cdot 10^{%$(round(Int, log10(enhance)))}$", fontsize = 15, padding = (150, 0, -300, 0))
    end
    Colorbar(fig[1:2, 3], colormap = :balance, limits = (-maxG, maxG), ticks = ([-maxG, maxG], niceticklabel.([-maxG, maxG])), label = labels["Filter"].barlabel, labelpadding = -25)
    colgap!(fig.layout, 1, 10)
    colgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 1, 10)

    Label(fig[1, 1:2, Top()], L"$L = %$(system.chain_params.L |> Int)$ nm, $\tau = %$(τ)$ ", padding = (0, 0, 5, 0))
    return fig
end

fig = fig_GVz("F_N=100"; maxG = 1e-2, enhance = 1e3)

##
Ns = [10, 50, 100, 500, 1000]
maxGs = [1e-1, 1e-1, 1e-2, 1e-2, 1e-2]
enhances = [1, 1e2, 1e3, 1e3, 1e3]
for (N, maxG, enhance) in zip(Ns, maxGs, enhances)
    fig = fig_GVz("F_N=$(N)"; maxG, enhance)
    save("plots/figures/F_N=$(N).pdf", fig)
end

##
N = 100
τ = 0.1
fig = fig_GVz("F_N=$(N)"; maxG = 1e-4, enhance = 1e5, τ)
save("plots/figures/F_N=$(N)_τ=$(τ).pdf", fig)
fig

## IRL - ILR
function I(G, dω)
    return sum(G, dims = 2) * dω
end
function fig_difIs(name::String)
    res = load("data/Conductance/$(name).jld2")["res"];
    @unpack system, Gs, path = res;
    @unpack params = system;
    @unpack ωrng, Vzrng, τrng = params;

    γ = system.NH_params.γdict[1][2]
    dω = real(ωrng[2] - ωrng[1])
    τrng = [0.1, 1]
    is = [1, 10]
    fig = Figure()
    xticks = ([first(Vzrng), γ, last(Vzrng)], [L"0", L"\gamma_y", L"%$(last(Vzrng))"])
    ax = Axis(fig[1, 1]; xlabel = L"$V_Z$ (meV)", ylabel = L"(I_{LR} - I_{RL})/I_{LR}", xticks)
    linestyle = :solid
    linewidth = 4
    color = :navyblue
    for (i, τ) in zip(is, τrng)
        IRL = sum(Gs[1, 2][:, :, i], dims = 1) * dω |> vec
        ILR = sum(Gs[2, 1][:, :, i], dims = 1) * dω |> vec
        ΔI = (IRL .- ILR) ./ IRL
        if i != 1
            linestyle = :dash
            linewidth = 2
            color = :red
        end
        lines!(ax, Vzrng, ΔI; linestyle, linewidth, color, label = L"$\tau = %$(τ)$")
    end
    axislegend(ax, position = :rb, frame = false)
    xlims!(ax, (first(Vzrng), last(Vzrng)))
    return fig
end
fig = fig_difIs("F_N=100")
save("plots/figures/F_N=100_difIs.pdf", fig)
fig

## Isolate state channel
function find_Eindex(Es, ωrng)
    return map(E -> findmin(abs.(E .- ωrng))[2], Es)
end

function find_Δ(G, iE, Vzrng, iτ; ia = 1, ib = 2, dωi::Int = 0)
    if dω == 0
    Δ = [
        abs((G[ia, ib][iE[i], i, iτ] - G[ib, ia][iE[i], i, iτ])/G[ia, ib][iE[i], i, iτ]) for i in eachindex(Vzrng)
    ]
    else 
        Δ = map(eachindex(Vzrng)) do i
            irng = (iE[i] - dωi):(iE[i] + dωi)
            return sum((G[ia, ib][irng, i, iτ] .- G[ib, ia][irng, i, iτ])) / sum(G[ia, ib][irng, i, iτ])
        end
    end
    return Δ
end
function fig_difGs(name::String; τ = 1)
    res = load("data/Conductance/$(name).jld2")["res"];
    eres = load("data/Spectrum/$(name).jld2")["res"];
    @unpack Es = eres;
    @unpack system, Gs, path = res;
    @unpack params = system;
    @unpack ωrng, Vzrng, τrng = params;

    ωrng = ωrng .|> real
    iτ = findmin(abs.(τrng .- τ))[2]
    dω0 = ωrng[2] - ωrng[1]
    γ = system.NH_params.γdict[1][2]

    dω = γ/2
    dωi = round(Int, dω/dω0)

    rEs = Es[:, iτ] .|> real
    rEs = [sort(e) for e in rEs]

    i1 = find_Eindex(getindex.(rEs, 1), ωrng)
    i2 = find_Eindex(getindex.(rEs, 2), ωrng)

    Δ1 = find_Δ(Gs, i1, Vzrng, iτ; dωi)
    Δ2 = find_Δ(Gs, i2, Vzrng, iτ; dωi)

    Δlocal = find_Δ(Gs, i1, Vzrng, iτ; ia = 1, ib = 1)

    fig = Figure()
    xticks = ([first(Vzrng), γ, last(Vzrng)], [L"0", L"\gamma_y", L"%$(last(Vzrng))"])

    ax = Axis(fig[1, 1]; xlabel = L"$V_Z$ (meV)", ylabel = L"$\left|(G_{LR} - G_{RL})/G_{LR} \right|$", xticks, )
    lines!(ax, Vzrng, Δ1; color = :purple, label = L"$1$", linewidth = 4)
    lines!(ax, Vzrng, Δ2; color = :purple, linestyle = :dash, label = L"$2$", linewidth = 4)
    lines!(ax, Vzrng, Δlocal; color = :red, linestyle = :solid, label = L"$(G_{RR} - G_{LL})/G_{RR}$", linewidth = 4)
    xlims!(ax, (first(Vzrng), last(Vzrng)))

    axislegend(ax, position = (1, 0.3),)
    return fig 
end

fig = fig_difGs("F_N=100"; τ = 1)
save("plots/figures/F_N=100_difGs.pdf", fig)
fig

## GLR GRL for a state
function get_GE(Gs, ia, ib, iE, iτ, Vzrng)
    G = map(eachindex(Vzrng)) do i
        return Gs[ia, ib][iE[i], i, iτ] |> abs
    end
    return G
end
function fig_GVzE(name::String; τ = 1)
    res = load("data/Conductance/$(name).jld2")["res"];
    eres = load("data/Spectrum/$(name).jld2")["res"];
    @unpack Es = eres;
    @unpack system, Gs, path = res;
    @unpack params = system;
    @unpack ωrng, Vzrng, τrng = params;

    ωrng = ωrng .|> real
    iτ = findmin(abs.(τrng .- τ))[2]

    rEs = Es[:, iτ] .|> real
    rEs = [sort(e) for e in rEs]

    i1 = find_Eindex(getindex.(rEs, 1), ωrng)
    i2 = find_Eindex(getindex.(rEs, 2), ωrng)

    fig = Figure()
    xticks = ([first(Vzrng), γ, last(Vzrng)], [L"0", L"\gamma_y", L"%$(last(Vzrng))"])
    ax = Axis(fig[1, 1]; xlabel = L"$V_Z$ (meV)", ylabel = L"$|G|$", xticks, yscale = log10 )

    G12_1 = get_GE(Gs, 1, 2, i1, iτ, Vzrng)
    G12_2 = get_GE(Gs, 1, 2, i2, iτ, Vzrng)
    G21_1 = get_GE(Gs, 2, 1, i1, iτ, Vzrng)
    G21_2 = get_GE(Gs, 2, 1, i2, iτ, Vzrng)

    
    lines!(ax, Vzrng, G12_1; color = :darkgreen, label = L"$G_{LR}^{(1)}$", linewidth = 4)
    lines!(ax, Vzrng, G12_2; color = :darkgreen, linestyle = :dash, label = L"$G_{LR}^{(2)}$", linewidth = 4)
    lines!(ax, Vzrng, G21_1; color = :navyblue, label = L"$G_{RL}^{(1)}$", linewidth = 4)
    lines!(ax, Vzrng, G21_2; color = :navyblue, linestyle = :dash, label = L"$G_{RL}^{(2)}$", linewidth = 4)

    xlims!(ax, (first(Vzrng), last(Vzrng)))
    axislegend(ax, position = (1, 1),)
    return fig 
end

fig = fig_GVzE("F_N=100"; τ = 1)
save("plots/figures/F_N=100_GVzE.pdf", fig)
fig