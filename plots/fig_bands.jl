function fig_bands(name::String; Vz = 0.2,  dir = "data/Conductance", colormap = :RdBu_9)
    path = "$(dir)/$(name).jld2"
    res = load(path)["res"]
    @unpack t, α, μ, a0 = res.system.chain_params
    γ0, γy = res.system.NH_params.γdict |> first |> last
    krng = range(-3π/4, 3π/4, length = 1000)
    Erng = range(-.7, 1.2, length = 1000)
    E(k, p, VZ) = γ0*1im - µ - cos(k) + 1 - p * sqrt(complex(VZ^2 + complex(γy*1im + (α/2a0) * sin(k))^2))

    fig = Figure(size = (600, 250), fontsize = 16)
    ax = Axis(fig[1, 1], 
        xlabel = L"k/a_0", 
        ylabel = L"E", 
        xticks = ([first(krng), 0, last(krng)],  [L"-\pi\leftarrow ", L"0", L"\rightarrow \pi"]),
        xticksize = 0,
        yticks = ([0], [L"\mu"])
    )

    function Ecolor(k, ω; p, Vz, Δω = 0.05)
        Eplus = E(k, p, Vz)
        if abs(ω - real(Eplus)) < Δω
            return imag(Eplus)
        else
            return NaN
        end
    end

    pts = Iterators.product(krng, Erng) 
    imap = map(pts) do (k, ω)
        Ecolor(k, ω; p = 1, Vz)
    end

    imap2 = map(pts) do (k, ω)
        Ecolor(k, ω; p = -1, Vz)
    end

    hmap = heatmap!(ax, krng, Erng, imap; colormap, )
    hmap = heatmap!(ax, krng, Erng, imap2; colormap, )
    #lines!(ax, krng, E.(krng, 1, Vz) .|> real; color = :black, linewidth = 1.5)
    #lines!(ax, krng, E.(krng, -1, Vz) .|> real; color = :black, linewidth = 1.5)
    #hlines!(ax, 0; color = :black, linestyle = :dash, linewidth = 1.5)
    ylims!(ax, first(Erng), last(Erng))

    kc = krng[findmin(abs.(real.(E.(krng, 1, Vz))))[2]]

    arrows!(ax, [kc - 0.1], [-0.15], [-.27], [0.3]; linewidth = 3, color = cgrad(colormap)|>last)
    arrows!(ax, [-kc + 0.1], [-0.15], [.27], [0.3]; linewidth = 3, color = cgrad(colormap)|>first)


    Colorbar(fig[1, 2]; colormap, label = L"\text{Im}(E)", ticks = ([0, 1], [L"\rightarrow 0", L"\gamma"]), labelpadding = -30)
    colgap!(fig.layout, 1, 5)

    Label(fig[1, 1, TopLeft()], "a"; padding = (-40, 0, -25, 0), labstyle...)

    return fig
end

fig = fig_bands("base_wire_nh")
save("plots/figures/fig_bands.pdf", fig)
fig