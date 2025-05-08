using Pkg
Pkg.activate("plots")
Pkg.instantiate()
using CairoMakie, JLD2, Parameters

## Conductance figures

function plot_conductance(pos, Gs, ωrng, μrng; colorrange = (-.1, .1))
    ax = Axis(pos; xlabel = L"\mu / t", ylabel = L"\omega / t")
    hmap = heatmap!(ax, μrng, real.(ωrng), Gs; colormap = :balance, colorrange)
    #vlines!(ax, 1; color = :darkgreen, linestyle = :dash)
    return ax, hmap
end

function plot_over_spectrum(ax, µrng, Es)
    for E in eachrow(Es)
        #scatter!(ax, µrng[1:5:end], real.(E[1:5:end]); color = (:white, 0.5), markersize = 2)
        scatter!(ax, µrng, real.(E); color = (:green, 0.5), markersize = 1)
        scatter!(ax, µrng, imag.(E); color = (:red, 0.5), markersize = 1)
    end
end

contact_dict = Dict(
    1 => "R",
    2 => "L",
)

function fig_conductance(name::String; maxG = 0.05, trans_coef = 0.01)
    res = load("data/Conductance/$(name).jld2")["res"]
    eres = load("data/Spectrum/$(name).jld2")["res"]
    @unpack Es = eres
    @unpack system, Gs, path = res
    @unpack params = system
    @unpack ωrng, µrng = params

    fig = Figure()
    
    for i in axes(Gs, 2), j in axes(Gs, 1)
        lim = maxG
        if i != j
            lim = trans_coef * maxG
        end
        ax, hmap = plot_conductance(fig[i, j], Gs[i, j]', ωrng, μrng; colorrange = (-lim, lim))
        i == j && plot_over_spectrum(ax, eres.system.params.μrng, Es)
        vlines!(ax, 1; color = :darkgreen, linestyle = :dash)
        ylims!(ax, (-2.5, 2.5))
        xlims!(ax, (0, 4))

        j == 2 && hideydecorations!(ax, ticks = false, grid = false)
        i == 1 && hidexdecorations!(ax, ticks = false, grid = false)

        Label(fig[i, j, Top()], L"$G_{%$(contact_dict[i]) %$(contact_dict[j])}$", fontsize = 15, padding = (180, 0, -140, 0))

    end
    Colorbar(fig[1:2, 3], colormap = :balance, limits = (-maxG, maxG), ticks = [-maxG, maxG], label = L"G", labelpadding = -25)
    colgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 1, 5)
    return fig
end

fig = fig_conductance("hermitian_chain")
save("plots/figures/conductance_hermitian_chain.pdf", fig)
fig

##
fig = fig_conductance("nh_even"; trans_coef = 0.0001)
save("plots/figures/conductance_nh_even.pdf", fig)
fig

##
fig = fig_conductance("nh_odd"; trans_coef = 1, maxG = 1)
save("plots/figures/conductance_nh_odd.pdf", fig)
fig

##
fig = fig_conductance("nh_odd_left"; trans_coef = 0.01)
save("plots/figures/conductance_nh_odd_left.pdf", fig)
fig