# Period fit with the period as a free parameter (FALC / CdR style), Fourier-order scan, phase ephemeris.

"""
    PeriodFit

Result of [`fit_period`](@ref) (or [`lightcurve_fit`](@ref) at a fixed period). Fields: `period`, `period_err` [same unit
as the times; NaN if the period was fixed], `order`, `a`, `b` (Fourier coefficients), `zero_points`, `groups` (group
labels in order), `labels` (display labels of the groups), `chi2`, `dof`, `chi2_red`, `amplitude` (peak-to-peak of the
Fourier series), `t0` (epoch of phase 0), `tmin` and `tmax` (first times after `t0` of minimum magnitude, i.e.
maximum brightness, and of maximum magnitude), and the data `t`, `y`, `σ`, `g` (group index per point). Plot it with `plot(fit)`
(phased light curve in the style of the Minor Planet Bulletin; needs Plots.jl).
"""
struct PeriodFit
    period::Float64
    period_err::Float64
    order::Int
    a::Vector{Float64}
    b::Vector{Float64}
    zero_points::Vector{Float64}
    groups::Vector{Any}
    labels::Vector{String}
    chi2::Float64
    dof::Int
    chi2_red::Float64
    amplitude::Float64
    t0::Float64
    tmin::Float64
    tmax::Float64
    t::Vector{Float64}
    y::Vector{Float64}
    σ::Vector{Float64}
    g::Vector{Int}
    object::String
end

"""
    fourier_series(fit::PeriodFit, phase) -> Δmag

The fitted Fourier series (without zero point) at the rotational phase(s) `phase` (0–1, zero at `fit.t0`).
"""
fourier_series(fit::PeriodFit, ph::Real) = sum(fit.a[k] * cos(2π * k * ph) + fit.b[k] * sin(2π * k * ph) for k in 1:fit.order)
fourier_series(fit::PeriodFit, ph::AbstractVector) = fourier_series.(Ref(fit), ph)

"""
    phase_ephemeris(fit::PeriodFit, t) -> (phase, Δmag)

Rotational phase and model magnitude (Fourier series, without zero point) at the times `t`.
"""
function phase_ephemeris(fit::PeriodFit, t)
    ph = mod.((t .- fit.t0) ./ fit.period, 1.0)
    return ph, fourier_series(fit, ph)
end

"""
    residuals(fit::PeriodFit) -> Vector

Observed minus model magnitudes (zero points removed).
"""
function residuals(fit::PeriodFit)
    ph = mod.((fit.t .- fit.t0) ./ fit.period, 1.0)
    return fit.y .- fit.zero_points[fit.g] .- fourier_series(fit, ph)
end

function _build_fit(t, y, σ, P, Perr, order, groups, labels, t0, object)
    n = length(t)
    gidx, ng, u = _groups(groups, n)
    ff = falc_fit(t .- t0, y, 1 / P; σ=σ, order=order, groups=groups)
    dof = n - (ng + 2order) - (isnan(Perr) ? 0 : 1)
    phs = range(0, 1; length=20001)
    m = [sum(ff.a[k] * cos(2π * k * p) + ff.b[k] * sin(2π * k * p) for k in 1:order) for p in phs]
    φ_magmin = phs[argmin(m)]                                    # minimum magnitude = maximum brightness
    φ_magmax = phs[argmax(m)]                                    # maximum magnitude = minimum brightness
    lab = labels === nothing ? (groups === nothing ? ["all"] : string.(u)) : collect(String.(labels))
    length(lab) == ng || throw(DimensionMismatch("labels must have one entry per group"))
    return PeriodFit(P, Perr, order, ff.a, ff.b, ff.zero_points, collect(u), lab, ff.chi2, dof, ff.chi2 / max(dof, 1),
                     ff.amplitude, t0, t0 + P * φ_magmin, t0 + P * φ_magmax, collect(Float64, t), collect(Float64, y),
                     collect(Float64, σ === nothing ? ones(n) : σ), gidx, object)
end

"""
    lightcurve_fit(t, y, P; σ=nothing, order=4, groups=nothing, labels=nothing, t0=minimum(t), object="") -> PeriodFit

Fourier fit of the light curve at the fixed period `P` (no period error), e.g. to plot a light curve folded at a
trial period.
"""
lightcurve_fit(t, y, P; σ=nothing, order::Integer=4, groups=nothing, labels=nothing, t0=minimum(t), object="") =
    _build_fit(t, y, σ, P, NaN, order, groups, labels, t0, object)

"""
    fit_period(t, y, P0; σ=nothing, order=4, groups=nothing, labels=nothing, t0=minimum(t), search=0.02, object="") -> PeriodFit

Least-squares fit of a Fourier series of the given `order` with the period free, starting from `P0` (same unit as
`t`; for asteroid light curves `t` in days, light-time corrected) and one zero point per group (session, night, set of
comparison stars). The Fourier coefficients and zero points are solved linearly at each trial period (variable
projection); the period minimises χ² within `P0 (1 ± search)`. The period error follows from the curvature of χ²(P),
scaled by the reduced χ² (as for an unweighted or poorly weighted fit; Harris et al. 1989). Equivalent to the CdR / FALC
fit with the period as a free parameter.
"""
function fit_period(t, y, P0::Real; σ=nothing, order::Integer=4, groups=nothing, labels=nothing, t0=minimum(t),
                    search::Real=0.02, object="")
    χ2(P) = first(falc(t .- t0, y, [1 / P]; σ=σ, order=order, groups=groups))[1]
    # golden-section search on a bracketing grid
    Ps = range(P0 * (1 - search), P0 * (1 + search); length=401)
    c = [χ2(P) for P in Ps]
    k = argmin(c)
    lo, hi = Ps[max(k - 1, 1)], Ps[min(k + 1, end)]
    gr = (sqrt(5) - 1) / 2
    x1 = hi - gr * (hi - lo); x2 = lo + gr * (hi - lo); f1 = χ2(x1); f2 = χ2(x2)
    for _ in 1:80
        if f1 < f2
            hi = x2; x2 = x1; f2 = f1; x1 = hi - gr * (hi - lo); f1 = χ2(x1)
        else
            lo = x1; x1 = x2; f1 = f2; x2 = lo + gr * (hi - lo); f2 = χ2(x2)
        end
        hi - lo < 1e-12 * P0 && break
    end
    P = (lo + hi) / 2
    T = maximum(t) - minimum(t)
    h = max(1e-4 * P^2 / T, 1e-9 * P)
    d2 = (χ2(P + h) - 2χ2(P) + χ2(P - h)) / h^2
    n = length(t); gidx, ng, _ = _groups(groups, n)
    χ2r = χ2(P) / max(n - ng - 2order - 1, 1)
    Perr = d2 > 0 ? sqrt(2 / d2 * χ2r) : NaN                    # covariance scaled by χ²_red, as scipy curve_fit
    (k == 1 || k == length(Ps)) && @warn "fit_period: χ² minimum at the edge of the search interval; widen `search`"
    return _build_fit(t, y, σ, P, Perr, order, groups, labels, t0, object)
end

"""
    fourier_order_scan(t, y, P; σ=nothing, groups=nothing, orders=1:10) -> Vector{NamedTuple}

χ², degrees of freedom, reduced χ², amplitude and BIC of Fourier fits of increasing `order` at the period `P` (Harris
et al. 1989 recommend the lowest order that fits the data; the minimum BIC is a guide).
"""
function fourier_order_scan(t, y, P; σ=nothing, groups=nothing, orders=1:10)
    n = length(t); _, ng, _ = _groups(groups, n)
    map(collect(orders)) do m
        f = falc_fit(t, y, 1 / P; σ=σ, order=m, groups=groups)
        k = ng + 2m; dof = n - k
        (order=m, chi2=f.chi2, dof=dof, chi2_red=f.chi2 / max(dof, 1), amplitude=f.amplitude, bic=f.chi2 + k * log(n))
    end
end

function Base.show(io::IO, f::PeriodFit)
    @printf(io, "PeriodFit(%sP = %.6f ± %.6f (%.5f ± %.5f h), order %d, amplitude %.3f, χ²_red %.2f, %d points, %d groups)",
            isempty(f.object) ? "" : f.object * ": ", f.period, f.period_err, 24f.period, 24f.period_err, f.order,
            f.amplitude, f.chi2_red, length(f.t), length(f.groups))
end
