function sumψ(ψ)
    ψodd = ψ[1:2:end]
    ψeven = ψ[2:2:end]
    return 0.5*(abs.(ψodd).^2 + abs.(ψeven).^2) |> vec
end

function separate(Ψ, n)
    ΨRo = Ψs[1][:, n][1:2:end]
    ΨRe = Ψs[1][:, n][2:2:end]
    ΨLo = Ψs[2][:, n][1:2:end]
    ΨLe = Ψs[2][:, n][2:2:end]
    return ΨRo, ΨRe, ΨLo, ΨLe
end

function project(ψodd, ψeven, v)
    return [dot(v, [Ψo, Ψe]) for (Ψo, Ψe) in zip(ψodd, ψeven)]
end

function fig_bands(name::String; nev = 2, Vz = 0.2,  dir = "data/Conductance", dirS = "data/States", colormap = cgrad(:PuOr_4) |> reverse, colorpsi = cgrad(:managua10)[2], colorL = cgrad(:managua10)[8], )
    path = "$(dir)/$(name).jld2"
    res = load(path)["res"]
    @unpack t, α, μ, a0 = res.system.chain_params

    pathS = "$(dirS)/base_wire_nh_long.jld2"
    resS = load(pathS)["res"]
    @unpack Ψs = resS

    γ0, γy = res.system.NH_params.γdict |> first |> last
    krng = range(-3π/4, 3π/4, length = 1000)
    Erng = range(-.7, 1.2, length = 1000)
    E(k, p, VZ) = γ0*1im - µ - cos(k) + 1 - p * sqrt(complex(VZ^2 + complex(γy*1im + (α/2a0) * sin(k))^2))

    fig = Figure(size = (600, 500), fontsize = 16)
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

    colordown = cgrad(colormap)|> last 
    colorup = cgrad(colormap)|> first

    arrows!(ax, [kc - 0.1], [-0.15], [-.27], [0.3]; linewidth = 3, color = colordown)
    arrows!(ax, [-kc + 0.1], [-0.15], [.27], [0.3]; linewidth = 3, color = colorup)


    Colorbar(fig[1, 2]; colormap, label = L"\text{Im}(E)", ticks = ([0, 1], [L"\rightarrow 0", L"\gamma"]), labelpadding = -30)
    colgap!(fig.layout, 1, 5)


    ax = Axis(fig[2, 1]; xlabel = L"x", ylabel = L"$\left| \psi \right|^2$ (arb. units)")
    N = length(Ψs[1][:, 1])/2 
    ax.xticks = ([0, N], [L"0", L"L"])
    ax.yticks = ([0, 1], [L"0", L"1"])
    for n in 2:nev
        ΨRo, ΨRe, ΨLo, ΨLe = separate(Ψs, n)
        ΨRu, ΨRd = map(v -> project(ΨRo, ΨRe, v), [[1im, 1], [-1im, 1]] ./sqrt(2))
        ΨLu, ΨLd = map(v -> project(ΨLo, ΨLe, v), [[1im, 1], [-1im, 1]] ./sqrt(2))

        ΨRu2 = ΨRu .|> norm 
        ΨRd2 = ΨRd .|> norm
        nu = -ΨLu .* ΨRu .|> real
        nd = ΨLd .* ΨRd .|> real

        mΨs = maximum(vcat(ΨRu2, ΨRd2))
        mns = maximum(vcat(nu, nd))
        lines!(ax, ΨRu2 ./ mΨs; color = colorpsi)
        lines!(ax, ΨRd2 ./ mΨs; color = colorpsi, linestyle = :dash)
        #lines!(ax, ΨLu .|> norm; color = :red)
        #lines!(ax, ΨLd .|> norm; color = :red, linestyle = :dash)
        lines!(ax, nu ./ mns; color = colorL)
        lines!(ax, nd ./ mns; color = colorL, linestyle = :dash)
    end

    elem_1 = [LineElement(color = colorpsi, )]
    elem_2 = [LineElement(color = colorL, )]
    elem_3 = [LineElement(color = :black, linestyle = :solid)]
    elem_4 = [LineElement(color = :black, linestyle = :dash)]

    axislegend(ax, 
        [elem_1, elem_2, elem_3, elem_4], 
        [L"\left|\psi_1^r\right|^2", L"\text{LDOS}", L"\downarrow_y", L"\uparrow_y"],
        position = (0.5, 0.8),
        labelpadding = -20,
        framevisible = false,
    )

    xlims!(ax, (0, N))
    ylims!(ax, (0, 1.1))

    Label(fig[1, 1, TopLeft()], "a"; padding = (-40, 0, -25, 0), labstyle...)
    Label(fig[2, 1, TopLeft()], "b"; padding = (-40, 0, -25, 0), labstyle...)

    return fig
end

fig = fig_bands("base_wire_nh")
save("plots/figures/fig_bands.pdf", fig)
fig