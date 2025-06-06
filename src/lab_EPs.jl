using CairoMakie

γy = 1
α = 1

d(Vz, k) = sqrt(complex(Vz^2 + (α * k + 1im * γy)^2))

Vrng = range(0, 2, length = 201)
k = 0 + 0.5im

ds = d.(Vrng, k)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$V_z$", ylabel = L"$d$")
lines!(ax, Vrng, ds .|> real; color = :navyblue, linewidth = 2)
lines!(ax, Vrng, ds .|> imag; color = :red, linewidth = 2)
fig