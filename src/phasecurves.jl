# Phase curves: reduced magnitudes, linear phase-curve fits with band offsets, photometric transformations.

"""
    reduced_magnitude(m, r, Δ) -> m − 5 log10(r Δ)

Magnitude reduced to unit heliocentric (`r`) and geocentric (`Δ`) distance [au].
"""
reduced_magnitude(m, r, Δ) = m .- 5 .* log10.(r .* Δ)

"""
    fit_phase_curve(α, mred; σ=nothing, bands=nothing, ref=nothing, beta=nothing, clip=3.0, maxiter=10) -> NamedTuple

Linear phase curve `mred = H + β α + offset(band)` by weighted least squares with iterative `clip`·σ clipping (σ
rescaled by the reduced χ² when it exceeds 1). `α` is the solar phase angle [deg], `mred` the reduced magnitude
([`reduced_magnitude`](@ref)). With `bands` (labels per point) one offset per band relative to the reference band `ref`
(default: the first label) is fitted, so `H` refers to `ref`. `beta` fixes the phase coefficient [mag/deg].
Returns `(H, H_err, beta, beta_err, offsets::Dict, offset_errs::Dict, keep, chi2_red, residuals)`; errors are scaled by
√max(χ²_red, 1). Suited to the small phase angles of distant objects (α ≲ 2°), where the H–G system is not appropriate.
"""
function fit_phase_curve(α::AbstractVector{<:Real}, mred::AbstractVector{<:Real}; σ=nothing, bands=nothing, ref=nothing,
                         beta=nothing, clip::Real=3.0, maxiter::Integer=10)
    n = length(α); length(mred) == n || throw(DimensionMismatch("α and mred differ in length"))
    w = _weights(σ, n); s = 1 ./ sqrt.(w)
    labels = bands === nothing ? fill("all", n) : collect(bands)
    ref = ref === nothing ? labels[1] : ref
    ref in labels || throw(ArgumentError("reference band $ref not present"))
    others = filter(!=(ref), unique(labels))
    design(k) = hcat(ones(count(k)), (beta === nothing ? (α[k],) : ())..., [Float64.(labels[k] .== b) for b in others]...)
    yv(k) = mred[k] .- (beta === nothing ? 0.0 : beta .* α[k])
    keep = trues(n)
    local c, C, χ2r
    for _ in 1:maxiter
        c, C, χ2 = _wls(design(keep), yv(keep), w[keep])
        dof = count(keep) - length(c)
        dof > 0 || throw(ArgumentError("too few points"))
        χ2r = χ2 / dof
        r = (yv(trues(n)) .- design(trues(n)) * c) ./ s
        knew = abs.(r) .< clip * sqrt(max(χ2r, 1.0))
        # a band whose points are all clipped would make the design singular: keep its points
        for b in others
            any(knew .& (labels .== b)) || (knew[labels .== b] .= keep[labels .== b])
        end
        knew == keep && break
        keep = knew
    end
    e = sqrt.(diag(C) .* max(χ2r, 1.0))
    nb = beta === nothing ? 2 : 1
    offs = Dict(b => c[nb+i] for (i, b) in enumerate(others)); offe = Dict(b => e[nb+i] for (i, b) in enumerate(others))
    offs[ref] = 0.0; offe[ref] = 0.0
    res = yv(trues(n)) .- design(trues(n)) * c
    return (H=c[1], H_err=e[1], beta=beta === nothing ? c[2] : Float64(beta), beta_err=beta === nothing ? e[2] : 0.0,
            offsets=offs, offset_errs=offe, keep=keep, chi2_red=χ2r, residuals=res, ref=ref)
end

"""
    ps1_to_V(r, g_r; g_r_err=0.0) -> (V, V_err)

Johnson V from Pan-STARRS1 r and g − r: V = r + 0.006 + 0.474 (g − r) (Tonry et al. 2012, ApJ 750, 99, Table 6;
scatter 0.012 mag, included in `V_err` together with the propagated colour error).
"""
ps1_to_V(r, g_r; g_r_err=0.0) = (r + 0.006 + 0.474g_r, sqrt((0.474g_r_err)^2 + 0.012^2))
