labels = Dict(
    "Wire" => (; xlabel = L"$B$ (meV)", ylabel = L"$V$ (mV)", barlabel = L"$G$ ($e^2/h$)"),
)

function plot_conductance(pos, Gs, ωrng, μrng; colorrange = (-.1, .1), labels = labels["Wire"])
    ax = Axis(pos; xlabel = labels.xlabel, ylabel = labels.ylabel)
    hmap = heatmap!(ax, μrng, real.(ωrng), Gs'; colormap = :balance, colorrange, rasterize = 5)
    #vlines!(ax, 1; color = :darkgreen, linestyle = :dash)
    return ax, hmap
end

function plot_conductance(pos, Gs, ωrng, xrng, yrng, y; colorrange = (-.1, .1), labels = labels["Wire"])
    iy = findmin(x -> abs(x - y), yrng)[2]
    Gs = Gs[:, :, iy]
    return plot_conductance(pos, Gs, ωrng, xrng; colorrange = colorrange, labels = labels)
end

function plot_over_spectrum(ax, µrng, Es; im = true, kw...)
    for E in eachrow(Es)
        #scatter!(ax, µrng[1:5:end], real.(E[1:5:end]); color = (:white, 0.5), markersize = 2)
        scatter!(ax, µrng, real.(E); color = (:green, 0.9), kw...)
        im && scatter!(ax, µrng, imag.(E); color = (:red, 0.5), markersize = 1)
    end
end
 
function plot_over_spectrum(ax, xrng, yrng, y, Es; im = true)
    iy = findmin(x -> abs(x - y), yrng)[2]
    Es = Es[:, iy]
    Es = hcat(Es...)
    sort!(Es, dims = 1, by = x -> abs(real(x)))
    plot_over_spectrum(ax, xrng, Es; im)
end

contact_dict = Dict(
    1 => "L",
    2 => "R",
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

labstyle = (font = "CMU Serif Bold", fontsize   = 20)
