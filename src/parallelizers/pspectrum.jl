function pspectrum(hfunc, xrng; nev = 20, ncv = nev, which = :LM, sigma = 0.0, check = 1, tol = 1e-10, kw...)
    function f(x; nev, kw...)
        try
            λ, ψ = eigs(hfunc(x); nev, kw...)
        catch
            return fill(NaN, nev)
        end
        if length(λ) < nev
            λ = vcat(λ, fill(NaN, nev - length(λ)))
        end
        return λ
    end
    Es = pfunction(f, [xrng]; nev, ncv, which, sigma, check, tol, kw...)
    Es = reshape(Es, size(xrng)...)
    Es = hcat(Es...)
    sort!(Es, dims = 1, by = x -> abs(real(x)))
    return Es
end

function pspectrum(hfunc, xrng, yrng; nev = 20, ncv = nev, which = :LM, sigma = 0.0, check = 1, tol = 1e-10, kw...)
    function f(x, y; nev, kw...)
        try
            λ, ψ = eigs(hfunc(x, y); nev, kw...)
        catch
            return fill(NaN, nev)
        end
        if length(λ) < nev
            λ = vcat(λ, fill(NaN, nev - length(λ)))
        end
        return λ
    end
    Es = pfunction(f, [xrng, yrng]; nev, ncv, which, sigma, check, tol, kw...)
    Es = reshape(Es, length(xrng), length(yrng))
    Es = hcat(Es...)
    sort!(Es, dims = 1, by = x -> abs(real(x)))
    return Es
end