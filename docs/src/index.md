# Lightcurves.jl

Rotational light curves and phase curves of small Solar System bodies: asteroids, Centaurs and trans-Neptunian objects.

The package is being rebuilt (towards version 0.1.0, November 2026). The current core covers the analysis of sparse and dense
photometry:

| Task | Functions |
|---|---|
| period search | [`gls`](@ref), [`falc`](@ref), [`falc_fit`](@ref), [`frequency_grid`](@ref), [`spectral_window`](@ref), [`permutation_fap`](@ref) |
| folded light curves | [`fold`](@ref), [`fourier_fit`](@ref), [`split_halves`](@ref) |
| phase curves | [`reduced_magnitude`](@ref), [`fit_phase_curve`](@ref), [`ps1_to_V`](@ref) |
| physical quantities | [`geometric_albedo`](@ref), [`diameter`](@ref), [`ellipsoid_amplitude`](@ref) |
| statistics | [`chi2_cdf`](@ref) |

Only the Julia standard library is required. Phase dispersion minimisation (PDM) is provided by the separate package
AstroPDM.jl.

```julia
using Lightcurves
f = frequency_grid(t; fmin=24 / 100, fmax=24 / 2)       # t in days → periods 2–100 h
p = gls(t, mag, f; σ=err)
χ2, χ0 = falc(t, mag, f; σ=err, order=2, groups=night)
P = 1 / f[argmin(χ2)]                                    # days
fit = falc_fit(t, mag, 1 / P; σ=err, groups=night)
fit.amplitude
```

Planned (see `STATUS.md` in the repository): H–G, H–G₁G₂ and H–G₁₂ phase functions, ALCDEF input and output, a
multi-apparition FALC with phase-angle term, AoV/ANOVA, ephemerides for the observation geometry, synthetic light curves
of shape models, and light-curve inversion.
