# Period search for unevenly sampled photometry: generalised Lomb–Scargle, FALC, spectral window.
# Times are in days (any origin), frequencies in cycles per day.

"""
    frequency_grid(t; fmin, fmax, oversample=10) -> Vector{Float64}

Evenly spaced trial frequencies from `fmin` to `fmax` (cycles per unit of `t`) with step `1 / (oversample * T)`, where `T`
is the time span of `t`. For rotation periods in hours, e.g. `fmin = 24 / 100, fmax = 24 / 2` with `t` in days.
"""
function frequency_grid(t::AbstractVector{<:Real}; fmin::Real, fmax::Real, oversample::Real=10)
    T = maximum(t) - minimum(t)
    T > 0 || throw(ArgumentError("the time span must be positive"))
    0 <= fmin < fmax || throw(ArgumentError("need 0 ≤ fmin < fmax"))
    return collect(fmin:1 / (oversample * T):fmax)
end

"""
    gls(t, y, f; σ=nothing) -> Vector{Float64}

Generalised Lomb–Scargle periodogram with a floating mean (Zechmeister & Kürster 2009, A&A 496, 577): normalised power
p(f) ∈ [0, 1] = 1 − χ²(f)/χ²₀ of a weighted sinusoid plus constant, for the trial frequencies `f`. `σ` are the
measurement errors (equal weights if omitted).
"""
function gls(t::AbstractVector{<:Real}, y::AbstractVector{<:Real}, f::AbstractVector{<:Real}; σ=nothing)
    n = length(t); length(y) == n || throw(DimensionMismatch("t and y differ in length"))
    w = _weights(σ, n); w = w ./ sum(w)
    ȳ = sum(w .* y); yy = sum(w .* (y .- ȳ) .^ 2)
    yy > 0 || throw(ArgumentError("y is constant"))
    P = similar(f, Float64)
    c = similar(t, Float64); s = similar(t, Float64)
    @inbounds for (j, fj) in enumerate(f)
        ω = 2π * fj
        C = S = YC = YS = CC = SS = CS = 0.0
        for i in 1:n
            ci, si = cos(ω * t[i]), sin(ω * t[i])
            wi = w[i]; dy = y[i] - ȳ
            C += wi * ci; S += wi * si; YC += wi * dy * ci; YS += wi * dy * si
            CC += wi * ci^2; SS += wi * si^2; CS += wi * ci * si
        end
        CC -= C^2; SS -= S^2; CS -= C * S
        D = CC * SS - CS^2
        P[j] = D > 0 ? (SS * YC^2 + CC * YS^2 - 2CS * YC * YS) / (D * yy) : 0.0
    end
    return P
end

# design matrix: one zero point per group + Fourier terms of order n at frequency f
function _falc_design(t, f, n, gidx, ng)
    m = length(t)
    X = zeros(m, ng + 2n)
    for i in 1:m
        X[i, gidx[i]] = 1.0
        ω = 2π * f * t[i]
        for k in 1:n
            X[i, ng+k] = cos(k * ω)
            X[i, ng+n+k] = sin(k * ω)
        end
    end
    return X
end
function _groups(groups, n)
    groups === nothing && return ones(Int, n), 1, [1]
    length(groups) == n || throw(DimensionMismatch("groups has the wrong length"))
    u = unique(groups)
    return [findfirst(==(g), u) for g in groups], length(u), u
end

"""
    falc(t, y, f; σ=nothing, order=2, groups=nothing) -> (χ², χ²₀)

FALC period scan (Harris et al. 1989, Icarus 77, 171): at each trial frequency a Fourier series of the given `order` is
fitted by weighted least squares together with one free zero point per group (`groups`: e.g. night, station or band
labels; one common zero point if omitted). Returns the χ² for each frequency and χ²₀ of the zero points alone; the
best frequency minimises χ². For a double-peaked asteroid light curve use `order ≥ 2` and expect the rotation
frequency, not twice it.
"""
function falc(t::AbstractVector{<:Real}, y::AbstractVector{<:Real}, f::AbstractVector{<:Real}; σ=nothing, order::Integer=2,
              groups=nothing)
    n = length(t); length(y) == n || throw(DimensionMismatch("t and y differ in length"))
    order >= 1 || throw(ArgumentError("order must be ≥ 1"))
    w = _weights(σ, n)
    gidx, ng, _ = _groups(groups, n)
    n > ng + 2order || throw(ArgumentError("too few points for $(ng) zero points and order $order"))
    X0 = zeros(n, ng); for i in 1:n; X0[i, gidx[i]] = 1.0; end
    _, _, χ0 = _wls(X0, y, w)
    χ2 = similar(f, Float64)
    for (j, fj) in enumerate(f)
        X = _falc_design(t, fj, order, gidx, ng)
        sw = sqrt.(w)
        A = sw .* X; b = sw .* y
        c = A \ b
        χ2[j] = sum(abs2, b .- A * c)
    end
    return χ2, χ0
end

"""
    falc_fit(t, y, f; σ=nothing, order=2, groups=nothing) -> NamedTuple

The FALC model at one frequency `f`: `(zero_points, groups, a, b, cov, chi2, amplitude, model)`, where `a[k]`, `b[k]` are
the cosine and sine coefficients, `amplitude` the peak-to-peak amplitude of the Fourier series and `model(phase)` the
series (without zero point) as a function of the rotational phase `phase = f t mod 1`.
"""
function falc_fit(t::AbstractVector{<:Real}, y::AbstractVector{<:Real}, f::Real; σ=nothing, order::Integer=2, groups=nothing)
    n = length(t)
    w = _weights(σ, n)
    gidx, ng, u = _groups(groups, n)
    X = _falc_design(t, f, order, gidx, ng)
    c, C, χ2 = _wls(X, y, w)
    a = c[ng+1:ng+order]; b = c[ng+order+1:end]
    model(ph) = sum(a[k] * cos(2π * k * ph) + b[k] * sin(2π * k * ph) for k in 1:order)
    m = model.(range(0, 1; length=1001))
    return (zero_points=c[1:ng], groups=u, a=a, b=b, cov=C, chi2=χ2, amplitude=maximum(m) - minimum(m), model=model)
end

"""
    spectral_window(t, f) -> Vector{Float64}

Spectral window |Σ exp(2πi f t)|² / N² of the sampling times `t` at the frequencies `f` (1 at f = 0). Peaks of a
periodogram that coincide with peaks of the window (e.g. 1 cycle per sidereal day, 1.0027 c/d) are aliases.
"""
spectral_window(t::AbstractVector{<:Real}, f::AbstractVector{<:Real}) =
    [abs2(sum(cis(2π * fj * ti) for ti in t)) / length(t)^2 for fj in f]

"""
    permutation_fap(stat, t, y, f; nperm=100, groups=nothing, rng=Random.default_rng(), kw...) -> Float64

False-alarm probability of the maximum of the periodogram function `stat(t, y, f; kw...)` (e.g. `gls`) by permutation:
`y` is shuffled `nperm` times (within each group, if `groups` is given) and the fraction of shuffles whose maximum
reaches the observed one is returned. Caution: shuffling single points ignores correlated errors (e.g. nightly zero
points); with sparse survey data the FAP is then too optimistic — shuffle nightly means or use `groups` = nights.
"""
function permutation_fap(stat, t, y, f; nperm::Integer=100, groups=nothing, rng=Random.default_rng(), kw...)
    obs = maximum(stat(t, y, f; kw...))
    gidx, ng, _ = _groups(groups, length(y))
    hits = 0
    ys = similar(y, Float64)
    for _ in 1:nperm
        copyto!(ys, y)
        for g in 1:ng
            idx = findall(==(g), gidx)
            ys[idx] = y[shuffle(rng, idx)]
        end
        hits += maximum(stat(t, ys, f; kw...)) >= obs
    end
    return hits / nperm
end
