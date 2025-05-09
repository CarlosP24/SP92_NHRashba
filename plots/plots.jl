using Pkg
Pkg.activate("plots")
Pkg.instantiate()
using CairoMakie, JLD2, Parameters, Quantica
include("../src/builders/Kitaev.jl")
include("../src/builders/Wire.jl")
## Conductance figures

function plot_conductance(pos, Gs, ωrng, μrng; colorrange = (-.1, .1), labels = labels["Kitaev"])
    ax = Axis(pos; xlabel = labels.xlabel, ylabel = labels.ylabel)
    hmap = heatmap!(ax, μrng, real.(ωrng), Gs; colormap = :balance, colorrange)
    #vlines!(ax, 1; color = :darkgreen, linestyle = :dash)
    return ax, hmap
end

function plot_over_spectrum(ax, µrng, Es; im = true)
    for E in eachrow(Es)
        #scatter!(ax, µrng[1:5:end], real.(E[1:5:end]); color = (:white, 0.5), markersize = 2)
        scatter!(ax, µrng, real.(E); color = (:green, 0.9), markersize = 1)
        im && scatter!(ax, µrng, imag.(E); color = (:red, 0.5), markersize = 1)
    end
end

contact_dict = Dict(
    1 => "R",
    2 => "L",
)

labels = Dict(
    "Kitaev" => (; xlabel = L"\mu / t", ylabel = L"\omega / t", barlabel = L"$G$ ($e^2/h$)"),
    "Wire_µ" => (; xlabel = L"$\mu$ (meV)", ylabel = L"$\omega$ (meV)", barlabel = L"$G$ ($e^2/h$)"),
    "Wire_Vz" => (; xlabel = L"$V_z$ (meV)", ylabel = L"$\omega$ (meV)", barlabel = L"$G$ ($e^2/h$)"),
)

function niceticklabel(num)
    anum = abs(num)
    if anum >= 0.01
        return L"%$(num)"
    end
    exp = round(Int, log10(anum))
    coef = ceil(Int, anum * 10^float(exp))
    return L"%$(sign(num) * coef |> Int) \cdot 10^{%$(exp)}"
end

function fig_conductance(name::String; maxG = 0.05, trans_coef = 0.01, ωlims = (-2.5, 2.5), im = true)
    res = load("data/Conductance/$(name).jld2")["res"]
    eres = load("data/Spectrum/$(name).jld2")["res"]
    @unpack Es = eres
    @unpack system, Gs, path = res
    @unpack params = system
    @unpack ωrng, x = params

    if x == :µ
        xrng = params.µrng
    elseif x == :Vz
        xrng = params.Vzrng
    else
        throw(ArgumentError("x must be :µ or :Vz"))
    end

    labs = labels
    if system.chain_params isa Wire_Params
        if x == :Vz
            labs = labels["Wire_Vz"]
        else
            labs = labels["Wire_µ"]
        end
    end

    fig = Figure()
    
    for i in axes(Gs, 2), j in axes(Gs, 1)
        lim = maxG
        if i != j
            lim = trans_coef * maxG
        end
        ax, hmap = plot_conductance(fig[i, j], Gs[i, j]', ωrng, xrng; colorrange = (-lim, lim), labels = labs)
        i == j && plot_over_spectrum(ax, xrng, Es; im)
        #vlines!(ax, 1; color = :darkgreen, linestyle = :dash)
        ylims!(ax, (first(ωrng) |> real, last(ωrng) |> real))
        xlims!(ax, (first(xrng), last(xrng)))

        j == 2 && hideydecorations!(ax, ticks = false, grid = false)
        i == 1 && hidexdecorations!(ax, ticks = false, grid = false)

        Label(fig[i, j, Top()], L"$G_{%$(contact_dict[i]) %$(contact_dict[j])}$", fontsize = 15, padding = (180, 0, -140, 0))

    end
    Colorbar(fig[1:2, 3], colormap = :balance, limits = (-maxG, maxG), ticks =( [-maxG, maxG], niceticklabel.([-maxG, maxG])), label = labs.barlabel, labelpadding = -25)
    colgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 1, 5)
    return fig
end

##
fig = fig_conductance("hermitian_chain")
save("plots/figures/conductance_hermitian_chain.pdf", fig)
fig

##
fig = fig_conductance("nh_even"; trans_coef = 0.0001)
save("plots/figures/conductance_nh_even.pdf", fig)
fig

##
fig = fig_conductance("nh_odd"; trans_coef = 1e4, maxG = .01)
save("plots/figures/conductance_nh_odd.pdf", fig)
fig

##
fig = fig_conductance("nh_odd_left"; trans_coef = 1e4, maxG = .01)
save("plots/figures/conductance_nh_odd_left.pdf", fig)
fig

##
fig = fig_conductance("nh_odd_right"; trans_coef = 1e4, maxG = .01)
save("plots/figures/conductance_nh_odd_right.pdf", fig)
fig

##
fig = fig_conductance("nh_odd_superleft"; trans_coef = 1e4, maxG = .01)
save("plots/figures/conductance_nh_odd_superleft.pdf", fig)
fig

##
fig = fig_conductance("Wire_base"; maxG = 1e-4, ωlims = (-0.25, 0.25), im = false, trans_coef = 1e-4 )
save("plots/figures/conductance_wire_base.pdf", fig)
fig

##
fig = fig_conductance("Wire_nh_101"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = false, trans_coef = 1e-6 )
save("plots/figures/conductance_wire_nh_101.pdf", fig)
fig

##
fig = fig_conductance("Wire_nh_001"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = false, trans_coef = 1e-4 )
save("plots/figures/conductance_wire_nh_001.pdf", fig)
fig

##
fig = fig_conductance("Wire_nh_100"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = false, trans_coef = 1e-8 )
save("plots/figures/conductance_wire_nh_100.pdf", fig)
fig