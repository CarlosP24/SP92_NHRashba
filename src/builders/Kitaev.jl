@with_kw struct Kitaev_Params
    μ::Float64 = 0.0
    t::Float64 = 1.0
    Δ::Float64 = 1.0 
    N::Int = 25
end

build(params::Kitaev_Params) = build_Kitaev(params)

function build_Kitaev(params::Kitaev_Params)
    @unpack µ, t, Δ, N = params

    lat = LP.linear(; ) |> supercell(N) |> supercell()

    kinetic = @onsite((; µ = 0) -> -µ * τz) + @hopping((; t = 0) -> -t * τz)
    pairing1 = @hopping((; Δ = 0) -> -Δ * 1im * τy;
    region = (r, dr) -> dr[1] > 0 )
    pairing2 = @hopping((; Δ = 0) -> Δ * 1im * τy;
    region = (r, dr) -> dr[1] < 0 )

    h = lat |> hamiltonian(kinetic + pairing1 + pairing2; orbitals = 2)
    return h 
end

@with_kw struct NH_Lead_Params
    Γodd::Float64 = 1.0
    Γeven::Float64 = 1.0
    Γodd_right = Γodd
    Γodd_left = Γodd
    Γeven_right = Γeven
    Γeven_left = Γeven
end

add_NH_lead!(params::NH_Lead_Params) = h -> add_NH_lead!(h, params)
function add_NH_lead!(h::Quantica.AbstractHamiltonian, params::NH_Lead_Params)
    @unpack Γodd_right, Γodd_left, Γeven_right, Γeven_left, dims = params

    nh_odd_right = @hopping!(
        (t, r, dr; Γodd_right = Γodd_right) -> t + Γodd_right * τ0;
        region = (r, dr) -> dr[1] > 0 )
    nh_odd_left = @hopping!(
        (t, r, dr; Γodd_left = Γodd_left) -> t - Γodd_left * τ0;
        region = (r, dr) -> dr[1] < 0 )

    nh_even_right = @hopping!(
        (t, r, dr; Γeven_right = Γeven_right) -> t + 1im*Γeven_right * τ0;
        region = (r, dr) -> dr[1] > 0 )
        
    nh_even_left = @hopping!(
        (t, r, dr; Γeven_left = Γeven_left) -> t + 1im*Γeven_left * τ0;
        region = (r, dr) -> dr[1] < 0 )

    if (Γodd_right == 0 && Γeven_right == 0) || (Γodd_left == 0 && Γeven_left == 0)
        return h
    elseif (Γodd_left == 0 && Γodd_right == 0)
        return h |> nh_even_left |> nh_even_right
    elseif (Γeven_left == 0 && Γeven_right == 0)
        return h |> nh_odd_right |> nh_odd_left
    else 
        return h |> nh_odd_right |> nh_odd_left |> nh_even_right |> nh_even_left
    end
end