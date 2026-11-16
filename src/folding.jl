# Folding, Fourier fits and the split-halves test.

"""
    fold(t, P; t0=minimum(t)) -> Vector{Float64}

Rotational phase in [0, 1) of the times `t` for the period `P` (same units), with zero at `t0`.
"""
fold(t::AbstractVector{<:Real}, P::Real; t0::Real=minimum(t)) =
    (P > 0 || throw(ArgumentError("P must be positive")); mod.((t .- t0) ./ P, 1.0))

"""
    fourier_fit(phase, y; σ=nothing, order=2) -> NamedTuple

Weighted least-squares Fourier series `c₀ + Σₖ aₖ cos(2πk φ) + bₖ sin(2πk φ)` in the rotational phase `φ`:
`(c0, a, b, cov, chi2, dchi2, amplitude, model)`, where `dchi2` is the χ² improvement over a constant and `model(φ)` the
fitted series including `c0`.
"""
function fourier_fit(phase::AbstractVector{<:Real}, y::AbstractVector{<:Real}; σ=nothing, order::Integer=2)
    n = length(phase); w = _weights(σ, n)
    X = hcat(ones(n), [cos.(2π * k .* phase) for k in 1:order]..., [sin.(2π * k .* phase) for k in 1:order]...)
    c, C, χ2 = _wls(X, y, w)
    χ0 = sum(w .* (y .- sum(w .* y) / sum(w)) .^ 2)
    a = c[2:order+1]; b = c[order+2:end]
    model(φ) = c[1] + sum(a[k] * cos(2π * k * φ) + b[k] * sin(2π * k * φ) for k in 1:order)
    m = model.(range(0, 1; length=1001))
    return (c0=c[1], a=a, b=b, cov=C, chi2=χ2, dchi2=χ0 - χ2, amplitude=maximum(m) - minimum(m), model=model)
end

"""
    split_halves(t, y, P; σ=nothing, nbins=10, t0=minimum(t), minpoints=2) -> NamedTuple

Split-halves test (Harris et al. 2014, Icarus 235, 55, Fig. 3): the light curve is folded at `P`, and the first
(φ < 0.5) and second half are binned in `nbins` bins of the phase within the half-cycle. Returns
`(chi2, dof, p, centers, mean1, err1, mean2, err2)`: χ² of the differences of the weighted bin means, its degrees of
freedom (bins with ≥ `minpoints` points in both halves) and the p-value. Similar halves (large p) mean that `P/2`
describes the data as well; different halves (small p) that the full period `P` is needed (e.g. two different maxima
of a double-peaked light curve). A symmetric double-peaked light curve gives similar halves too.
"""
function split_halves(t::AbstractVector{<:Real}, y::AbstractVector{<:Real}, P::Real; σ=nothing, nbins::Integer=10,
                      t0::Real=minimum(t), minpoints::Integer=2)
    n = length(t); w = _weights(σ, n)
    ph = fold(t, P; t0=t0)
    first_half = ph .< 0.5
    ψ = mod.(ph, 0.5) .* 2
    centers = [(b + 0.5) / nbins for b in 0:nbins-1]
    stats(k) = map(0:nbins-1) do b
        m = k .& (ψ .>= b / nbins) .& (ψ .< (b + 1) / nbins)
        count(m) < minpoints && return (NaN, NaN)
        ww = w[m]
        (sum(ww .* y[m]) / sum(ww), 1 / sqrt(sum(ww)))
    end
    s1 = stats(first_half); s2 = stats(.!first_half)
    m1, e1 = first.(s1), last.(s1); m2, e2 = first.(s2), last.(s2)
    ok = .!isnan.(m1) .& .!isnan.(m2)
    χ2 = sum(abs2, (m1[ok] .- m2[ok]) ./ sqrt.(e1[ok] .^ 2 .+ e2[ok] .^ 2))
    dof = count(ok)
    return (chi2=χ2, dof=dof, p=dof > 0 ? 1 - chi2_cdf(χ2, dof) : NaN, centers=centers, mean1=m1, err1=e1, mean2=m2, err2=e2)
end
