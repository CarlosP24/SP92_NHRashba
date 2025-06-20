function fig_spectrum(name::String; dirS = "data/Spectrum", dirG = "data/Conductance", maxG = 2e-2)
    path = "$(dir)/$(name).jld2"
    res = load(path)["res"]

    resG = load("$(dirG)/$(name).jld2")["res"]
    @unpack Gs, system = resG
    @unpack t, α, μ, a0 = res.system.chain_params
    t = 1
    α = 1
    a0 = 1
    @unpack Es, system = res
    @unpack Vzrng, ωrng = system.params
    ωrng = real.(ωrng)
    γ = system.NH_params.γdict |> first |> last |> first

    fig = Figure(size = (600, 500), fontsize = 16)

    ax, hmap = plot_conductance(fig[1, 1], Gs[1, 1], ωrng, Vzrng; colorrange = (-maxG, maxG))
    ax.yticks = ([-0.2, 0, 0.3], [L"-0.2", L"0", L"0.3"])
    ax.xticks = ([0, γ, 0.3], [L"0", L"γ_y", L"0.3"])
    ax.ylabelpadding = -15
    ax.xlabelpadding = -15
    for E in eachrow(Es)
        scatter!(ax, Vzrng, real.(E); color = (:green, 0.9), markersize = 3)
        #scatter!(ax, Vzrng, imag.(E); color = (:red, 0.5), markersize = 1)
    end
    ylims!(ax, (first(ωrng), last(ωrng)))
    vlines!(ax, γ; color = :black, linestyle = :dash)

    Colorbar(fig[1, 2]; colormap = cgrad(:balance, 256)[128:end], limits = (0, 1), label = L"$G_{\text{RR}}$ ($e^2/h$)", ticks = ([0, 1], [L"0", niceticklabel(maxG)]), labelpadding = -30)

    colgap!(fig.layout, 1, 5)

    fig_bands = fig[2, 1] = GridLayout()
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

    lines!(ax, krngp, E.(krngp, 1, Vz) .|> real; color = :red, linewidth = 1.5)
    lines!(ax, krngn, E.(krngn, -1, Vz) .|> real; color = :red, linewidth = 1.5)
    lines!(ax, krngp, E.(krngp, 1, Vz) .|> imag; color = :red, linestyle = :dash, linewidth = 1.5)
    lines!(ax, krngn, E.(krngn, -1, Vz) .|> imag; color = :red, linestyle = :dash, linewidth = 1.5)

    lines!(ax, krngp, E.(krngp, -1, Vz) .|> real; color = :blue, linewidth = 1.5)
    lines!(ax, krngn, E.(krngn, 1, Vz) .|> real; color = :blue, linewidth = 1.5)
    lines!(ax, krngp, E.(krngn, 1, Vz) .|> imag; color = :blue, linestyle = :dash, linewidth = 1.5)
    lines!(ax, krngn, E.(krngp, -1, Vz) .|> imag; color = :blue, linestyle = :dash, linewidth = 1.5)

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

    lines!(ax, krng, E.(krng, 1, Vz) .|> real; color = :red, linewidth = 1.5)
    lines!(ax, krng, E.(krng, -1, Vz) .|> real; color = :blue, linewidth = 1.5)

    lines!(ax, krng, E.(krng, 1, Vz) .|> imag; color = :red, linestyle = :dash, linewidth = 1.5)
    lines!(ax, krng, E.(krng, -1, Vz) .|> imag; color = :blue, linestyle = :dash, linewidth = 1.5)
    
    xlims!(ax, (first(krng), last(krng)))
    ylims!(ax, (-0.7, 2.0))

    hideydecorations!(ax, ticks = false, grid = false)

    colgap!(fig_bands, 1, 5)

    Label(fig_bands[1, 1, Top()], L"$B < \gamma_y$"; padding = (0, 0, 5, 0))
    Label(fig_bands[1, 2, Top()], L"$B > \gamma_y$"; padding = (0, 0, 5, 0))


    Label(fig[1, 1, TopLeft()], "a"; padding = (-40, 0, -25, 0), labstyle...)
    Label(fig_bands[1, 1, TopLeft()], "b"; padding = (-40, 0, 0, 0), labstyle...)
    Label(fig_bands[1, 2, TopLeft()], "c"; padding = (0, 0, 0, 0), labstyle...)

    rowgap!(fig.layout, 1, 5)

    return fig
end

fig = fig_spectrum("base_wire_nh")
save("plots/figures/fig_spectrum.pdf", fig)
fig