function kappa2(n, α, L, m, ħ2ome)
    θ = α * L * m / ħ2ome
    d = 1 - (θ / (n * π))^2
    f = sin(θ) / θ
    return (f / d)^2
end

function ep(n, Γ, params)
    @unpack α, L, m0, ħ2ome = params
    return Γ /  sqrt(kappa2(n, α, L, m0, ħ2ome))
end

function fig_spectrum(name::String; dirS = "data/Spectrum", dirG = "data/Conductance", maxG = 2e-2, colors = cgrad(:buda, 10, categorical = true)[[2, 8]])
    path = "$(dirS)/$(name).jld2"
    res = load(path)["res"]

    resG = load("$(dirG)/$(name).jld2")["res"]
    @unpack Gs, system = resG
    @unpack t, α, μ, a0 = system.chain_params
    t = 1
    α = 1
    a0 = 1
    @unpack Es, system = res
    Es = mapslices(x -> sort(x, by=real), Es; dims=1)
    @unpack Vzrng, ωrng = system.params
    ωrng = real.(ωrng)
    γ = system.NH_params.γdict |> first |> last |> first
    Vγ1 = ep(1, γ, system.chain_params)
    iVγ1 = findmin(abs.(Vzrng .- Vγ1))[2]

    fig = Figure(size = (600, 520), fontsize = 16)

    fig_cond = fig[2, 1] = GridLayout()

    ax, hmap = plot_conductance(fig_cond[1, 1], Gs[1, 1], ωrng, Vzrng; colorrange = (-maxG, maxG))
    ax.yticks = ([-0.2, 0, 0.3], [L"-0.2", L"0", L"0.3"])
    ax.xticks = ([0, 0.3], [L"0",  L"0.3"])
    ax.ylabelpadding = -15
    ax.xlabelpadding = -15
    for (i, E) in enumerate(eachrow(Es))
        color = get(colors, i, :green)
        if i != 2
            lines!(ax, Vzrng, real.(E); color, linewidth = 2)
        else
            lines!(ax, Vzrng[1:iVγ1], real.(E)[1:iVγ1]; color, linestyle = :dash)
            lines!(ax, Vzrng[iVγ1:end], real.(E)[iVγ1:end]; color)
        end
        #scatter!(ax, Vzrng, imag.(E); color = (colors[2], 0.5), markersize = 1)
    end
    ylims!(ax, (first(ωrng), last(ωrng)))
    vlines!(ax, γ; color = :black, linestyle = :dash, label = L"$\gamma_y$", linewidth = 2)

    vlines!(ax, Vγ1 ; ymin = 0.3, ymax = 0.6, color = :white, linestyle = :dash, label = L"$\gamma_y/\left|\kappa_1\right|$", linewidth = 2)

    vlines!(ax, [ep(2, γ, system.chain_params)]; ymin = 0.6, ymax = 0.9, color = :blue, linestyle = :dash, label = L"$\gamma_y/\left|\kappa_2\right|$", linewidth = 2)

    Colorbar(fig_cond[1, 2]; colormap = cgrad(:balance, 256)[128:end], limits = (0, 1), label = L"$G_{\text{RR}}$ ($e^2/h$)", ticks = ([0, 1], [L"0", niceticklabel(maxG)]), labelpadding = -30)

    axislegend( ax; position = (1, 0.6), orientation = :horizontal, padding = (0, 0, 10, 0), framevisible = false)
    colgap!(fig_cond, 1, 5)

    fig_bands = fig[1, 1] = GridLayout()
    krngp = range(0.001, 3π/4, length = 1000)
    krngn = range(-3π/4, -0.001, length = 1000)
    krng = range(-3π/4, 3π/4, length = 1000)
    γ = 1
    Vz = 0.9
    E(k, p, VZ) = γ*1im - µ - cos(k) + 1 - p * sqrt(complex(VZ^2 + complex(γ*1im + (α/2a0) * sin(k))^2))

    ax = Axis(fig_bands[1, 1],
        xlabel = L"k/a_0", 
        ylabel = L"E", 
        xticks = ([first(krng) + 0.5, 0, last(krng) - 0.5],  [L"-\pi\leftarrow ", L"0", L"\rightarrow \pi"]),
        xticksize = 0,
        yticks = ([0, γ], [L"\mu", L"\gamma_0"])
    )
    hidexdecorations!(ax, ticks = false, ticklabels = false, label = false)

    lines!(ax, krngp, E.(krngp, 1, Vz) .|> real; color = colors[1], linewidth = 1.5)
    lines!(ax, krngn, E.(krngn, -1, Vz) .|> real; color = colors[1], linewidth = 1.5)
    lines!(ax, krngp, E.(krngp, 1, Vz) .|> imag; color = colors[1], linestyle = :dash, linewidth = 1.5)
    lines!(ax, krngn, E.(krngn, -1, Vz) .|> imag; color = colors[1], linestyle = :dash, linewidth = 1.5)

    lines!(ax, krngp, E.(krngp, -1, Vz) .|> real; color = colors[2], linewidth = 1.5)
    lines!(ax, krngn, E.(krngn, 1, Vz) .|> real; color = colors[2], linewidth = 1.5)
    lines!(ax, krngp, E.(krngn, 1, Vz) .|> imag; color = colors[2], linestyle = :dash, linewidth = 1.5)
    lines!(ax, krngn, E.(krngp, -1, Vz) .|> imag; color = colors[2], linestyle = :dash, linewidth = 1.5)

    xlims!(ax, (first(krng), last(krng)))
    ylims!(ax, (-0.7, 2.0))

    axislegend(
        ax,
        [
            LineElement(color = :black, linestyle = :solid,),
            LineElement(color = :black, linestyle = :dash,)
        ],
        [L"\text{Re}", L"\text{Im}"],
        framevisible = false,
        orientation = :horizontal,
        position = (0.5, -.05)
    )


    Vz = 1.05
    ax = Axis(fig_bands[1, 2],
        xlabel = L"k/a_0", 
        ylabel = L"E", 
        xticks = ([first(krng) + 0.5, 0, last(krng) - 0.5],  [L"-\pi\leftarrow ", L"0", L"\rightarrow \pi"]),
        xticksize = 0,
        yticks = ([0, γ], [L"\mu", L"\gamma_0"])
    )

    hidexdecorations!(ax, ticks = false, ticklabels = false, label = false)

    lines!(ax, krng, E.(krng, 1, Vz) .|> real; color = colors[1], linewidth = 1.5)
    lines!(ax, krng, E.(krng, -1, Vz) .|> real; color = colors[2], linewidth = 1.5)

    lines!(ax, krng, E.(krng, 1, Vz) .|> imag; color = colors[1], linestyle = :dash, linewidth = 1.5)
    lines!(ax, krng, E.(krng, -1, Vz) .|> imag; color = colors[2], linestyle = :dash, linewidth = 1.5)
    
    xlims!(ax, (first(krng), last(krng)))
    ylims!(ax, (-0.7, 2.0))

    hideydecorations!(ax, ticks = false, grid = false)

    colgap!(fig_bands, 1, 5)

    Label(fig_bands[1, 1, Top()], L"$B < \gamma_y$"; padding = (0, 0, 5, 0))
    Label(fig_bands[1, 2, Top()], L"$B > \gamma_y$"; padding = (0, 0, 5, 0))


    Label(fig_bands[1, 1, TopLeft()], "a"; padding = (-40, 0, 0, 0), labstyle...)
    Label(fig_bands[1, 2, TopLeft()], "b"; padding = (0, 0, 0, 0), labstyle...)
    Label(fig_cond[1, 1, TopLeft()], "c"; padding = (-40, 0, -25, 0), labstyle...)

    rowgap!(fig.layout, 1, 5)

    return fig
end

fig = fig_spectrum("base_wire_nh")
save("plots/figures/fig_spectrum.pdf", fig)
fig