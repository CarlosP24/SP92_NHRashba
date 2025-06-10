function fig_conductance(name::String; y = 0.0, maxG = 0.05, trans_coef = 0.01, ωlims = (-2.5, 2.5), im = true, spectrum = true)
    res = load("data/Conductance/$(name).jld2")["res"]
    eres = load("data/Spectrum/$(name).jld2")["res"]
    @unpack Es = eres
    @unpack system, Gs, path = res
    @unpack params = system
    @unpack ωrng, x = params
    
    plot_c(pos, Gs, ωrng, xrng, yrng, y; kws...) = plot_conductance(pos, Gs, ωrng, xrng; kws...)
    plot_o(ax, xrng, yrng, y, Es; kws...) = plot_over_spectrum(ax, xrng, Es; kws...)
    yrng = nothing
    if x == :µ
        xrng = params.µrng
    elseif x == :Vz
        xrng = params.Vzrng
    elseif x == :θ
        xrng = params.θrng
    elseif x == (:Vz, :θ)
        xrng = params.Vzrng
        yrng = params.θrng
        plot_c = (pos, Gs, ωrng, xrng, yrng, y; kws...)  -> plot_conductance(pos, Gs, ωrng, xrng, yrng, y; kws...)
        plot_o = (ax, xrng, yrng, y, Es; kws...) -> plot_over_spectrum(ax, xrng, yrng, y, Es; kws...)
    else
        throw(ArgumentError("x must be :µ or :Vz"))
    end

    labs = labels["Kitaev"]
    if system.chain_params isa Wire_Params
        if x == :Vz
            labs = labels["Wire_Vz"]
        else
            labs = labels["Wire_µ"]
        end
    end

    if system.chain_params isa Filter_Params
        labs = labels["Filter"]
    end

    fig = Figure()
    
    for i in axes(Gs, 2), j in axes(Gs, 1)
        lim = maxG
        if i != j
            lim = trans_coef * maxG
        end
        ax, hmap = plot_c(fig[i, j], Gs[i, j], ωrng, xrng, yrng, y; colorrange = (-lim, lim), labels = labs)
        spectrum && i == j && plot_o(ax, xrng, yrng, y, Es; im)
        vlines!(ax, 0.05; color = :black, linestyle = :dash)
        #vlines!(ax, 1; color = :darkgreen, linestyle = :dash)
        ylims!(ax, (first(ωrng) |> real, last(ωrng) |> real))
        xlims!(ax, (first(xrng), last(xrng)))

        j == 2 && hideydecorations!(ax, ticks = false, grid = false)
        i == 1 && hidexdecorations!(ax, ticks = false, grid = false)

        Label(fig[i, j, Top()], L"$G_{%$(contact_dict[i]) %$(contact_dict[j])}$", fontsize = 15, padding = (-180, 0, -140, 0))

    end
    Colorbar(fig[1:2, 3], colormap = :balance, limits = (-maxG, maxG), ticks =( [-maxG, maxG], niceticklabel.([-maxG, maxG])), label = labs.barlabel, labelpadding = -25)
    
    #system.chain_params isa Filter_Params && Label(fig[1, 1:2, Top()], L"$\theta = %$(round(Int, y/π))\pi$")
    if occursin("=", name)
        τ_value = match(r"=(.*)", name)
        if τ_value !== nothing
            τ_num = tryparse(Float64, τ_value.captures[1])/10
            if τ_num !== nothing
            Label(fig[1, 1:2, Top()], L"\tau = %$(τ_num)")
            end
        end
    end
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
fig = fig_conductance("Wire_nh_101"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = true, trans_coef = 1e-6 )
save("plots/figures/conductance_wire_nh_101.pdf", fig)
fig

##
fig = fig_conductance("Wire_nh_001"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = true, trans_coef = 1e-4 )
save("plots/figures/conductance_wire_nh_001.pdf", fig)
fig

##
fig = fig_conductance("Wire_nh_100"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = true, trans_coef = 1e-8 )
save("plots/figures/conductance_wire_nh_100.pdf", fig)
fig

##
fig = fig_conductance("Wire_nh_110"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = true, trans_coef = 1e-8 )
save("plots/figures/conductance_wire_nh_110.pdf", fig)
fig

##
fig = fig_conductance("Wire_nh_010"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = true, trans_coef = 1e-4 )
save("plots/figures/conductance_wire_nh_010.pdf", fig)
fig

##
fig = fig_conductance("Wire_nh_001_strong"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = true, trans_coef = 1e-4 )
save("plots/figures/conductance_wire_nh_001.pdf", fig)
fig

##
fig = fig_conductance("Wire_nh_101_strong"; maxG = 5e-4, ωlims = (-0.25, 0.25), im = true, trans_coef = 1e-8 )
save("plots/figures/conductance_wire_nh_101.pdf", fig)
fig


##
fig = fig_conductance("Filter_base"; maxG = 1e-5, trans_coef = 1e-4)
fig

##
fig = fig_conductance("Filter_base_noSOC"; maxG = 1e-5,trans_coef = 1e-4)
fig

##
fig = fig_conductance("Filter_nh_example"; maxG = 1e-6, trans_coef = 1e-3)
fig

##
fig = fig_conductance("Filter_real"; maxG = 1e-5, trans_coef = 1e-5, spectrum = false)
fig

##
fig = fig_conductance("Filter_real_nh_0.05"; maxG = 1e-5, trans_coef = 1e-4, spectrum = false)
fig

##
fig = fig_conductance("Filter_real_nh_0.05_transparent"; maxG = 1, trans_coef = 1, spectrum = false)
fig

##
fig = fig_conductance("Filter_real_nh11_0.05"; maxG = 1e-5, trans_coef = 1e-6, spectrum = false)
fig

##
fig = fig_conductance("Filter_real_nh_0.01"; maxG = 1e-5, trans_coef = 1e-5, spectrum = false)
fig

##
fig = fig_conductance("Filter_real_nh11_0.01"; maxG = 1e-5, trans_coef = 1e-6, spectrum = false)
fig

##
fig = fig_conductance("Filter_real_nh_0.1"; maxG = 1e-7, trans_coef = 1e-5, spectrum = false)
fig

##
fig = fig_conductance("Filter_real_nh_0.5"; maxG = 1e-8, trans_coef = 1e-5)
fig

##
fig = fig_conductance("Filter_real_nh_1.0"; maxG = 1e-8, trans_coef = 1e-5)
fig

##
maxGs = [1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 1e-1, 1e-1, 1e-1, 1, 1]
for (τ, maxG) in zip(1:10, maxGs)
    fig = fig_conductance("Filter_real_nh_0.05_τ=$τ"; maxG, trans_coef = 1, spectrum = false)
    save("plots/figures/conductance_filter_real_nh_0.05_τ=$(τ)_clean.pdf", fig)
end

##
fig = fig_conductance("Filter_short_nh11_0.05_τ=0.1"; maxG = 1e-4, trans_coef = 1e-4, spectrum = true)