 labstyle = (font = "CMU Serif Bold", fontsize   = 20)
function fig_conductance(name::String; dir = "data/Conductance", maxGLL = 2e-2, maxGRL = 5e-5)
    path = "$(dir)/$(name).jld2"
    res = load(path)["res"]

    @unpack Gs, system = res
    @unpack Vzrng, ωrng = system.params

    γ = system.NH_params.γdict |> first |> last |> first
    ratio = maxGRL / maxGLL
    lratio = round(Int, log10(ratio))
    fig = Figure(size = (600, 500), fontsize = 16) 
    for i in 1:2, j in 1:2
        if i == j
            maxG = maxGLL
        else
            maxG = maxGRL
        end
        ax, hmap = plot_conductance(fig[i, j], Gs[i, j], ωrng, Vzrng; colorrange = (-maxG, maxG))
        vlines!(ax, γ; color = :black, linestyle = :dash)
        ax.xticks = ([0, γ, 0.2], [L"0", L"γ_y", L"0.2"])
        ax.yticks = ([-0.2, 0, 0.3], [L"-0.2", L"0", L"0.3"])
        ax.ylabelpadding = -15
        i == j && text!(ax, (0.15, -0.17); text = L"$G_{\text{%$(contact_dict[i])%$(contact_dict[j])}}$", align = (:center, :center), color = :black, fontsize = 16)
        i != j && text!(ax, (0.15, -0.17); text = L"$G_{\text{%$(contact_dict[i])%$(contact_dict[j])}} \cdot 10^{%$(-lratio)}$", align = (:center, :center), color = :black, fontsize = 16)
        i == 1 && hidexdecorations!(ax, ticks = false, grid = false)
        j != 1 && hideydecorations!(ax, ticks = false, grid = false)
        j != 1 && colgap!(fig.layout, j - 1, 15)
        i != 1 && rowgap!(fig.layout, i - 1, 5)
    end

    Colorbar(fig[1:2, 3], colormap = :balance, colorrange = (-maxGLL, maxGLL), label = L"$G$ ($e^2/h$)", ticks = ([-maxGLL, 0, maxGLL], [niceticklabel(-maxGLL), L"0", niceticklabel(maxGLL)]), labelpadding = -30)

    colgap!(fig.layout, 2, 5)

    Label(fig[1, 1, TopLeft()], "a"; padding = (-40, 0, -25, 0), labstyle...)
    Label(fig[1, 2, TopLeft()], "b"; padding = (-15, 0, -25, 0), labstyle...)
    Label(fig[2, 1, TopLeft()], "c"; padding = (-40, 0, -25, 0), labstyle...)
    Label(fig[2, 2, TopLeft()], "d"; padding = (-15, 0, -25, 0), labstyle...)
    return fig
end

fig = fig_conductance("base_wire_nh")
save("plots/figures/fig_conductance.pdf", fig)
fig