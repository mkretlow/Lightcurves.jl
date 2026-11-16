# Lightcurves.jl

Julia package for rotational light curves and phase curves of small Solar System bodies (asteroids, Centaurs,
trans-Neptunian objects). Part of a Julia toolkit for the study of small bodies (with Occultations.jl, AstroPDM.jl, …).

Version 0.1.0-DEV (rebuilt from scratch in November 2026; first release tag 0.1.0 planned before the light-curve inversion work). Light dependencies only (RecipesBase, JSON,
CodecZlib and standard libraries; FITSIO as optional extension for FITS images); plots through Plots.jl recipes.

## Features

| Task | Functions |
|---|---|
| Period search | `gls` (generalised Lomb–Scargle), `falc` / `falc_fit` (FALC Fourier-series scan with zero points per night/band), `frequency_grid`, `spectral_window`, `permutation_fap` |
| Folded light curves | `fold`, `fourier_fit`, `split_halves` (split-halves plot test) |
| Phase curves | `reduced_magnitude`, `fit_phase_curve` (linear phase curve with band offsets and σ clipping), `ps1_to_V` |
| Physical quantities | `geometric_albedo`, `diameter`, `ellipsoid_amplitude` |
| Period fit | `fit_period` (period free, error from χ² curvature, one zero point per session) → `PeriodFit`; `lightcurve_fit`, `fourier_order_scan`, `phase_ephemeris`, `residuals` |
| Files | `read_alcdef`, `write_alcdef`, `alcdef_table` (ALCDEF v2.3), `read_cdr` (CdR / CourbRot), `mpc_to_V` (MPC band corrections), `light_time_correct` |
| Data access | `fetch_alcdef_archive` / `read_alcdef_archive` / `alcdef_archive_objects` (local ALCDEF archive), `fetch_lcdb_archive` (local LCDB zip), `fetch_mpc_photometry` (MPC), `observing_geometry` / `fetch_horizons_geometry` (JPL Horizons), `fetch_sbdb` (JPL SBDB), `fetch_alcdef` / `alcdef_sessions` (ALCDEF), `lcdb_entry`, `fetch_lcdb_member`, `read_lcdb_summary` (LCDB) |
| CCD images | `masterframe`, `calibrate`, `findstars`, `centroid`, `aperture`, `instrumental_magnitude`, `fwhm_estimate`, `match_offset`, `platematch`, `robust_stats`; `readfits` with `using FITSIO` (package extension) |
| Plots (Plots.jl) | `plot(::PeriodFit)` (Minor Planet Bulletin style, optional residual panel), `plot(::Vector{ALCDEFBlock})`, `periodogramplot`, `phasecurveplot` (bands joint or per band), `splithalvesplot`, `windowplot` |

Phase dispersion minimisation (PDM): see the separate package AstroPDM.jl.

## Usage

```julia
using Pkg; Pkg.develop(path="/path/to/Lightcurves.jl")
using Lightcurves

f = frequency_grid(t; fmin=24 / 100, fmax=24 / 2)            # t in days → periods 2–100 h
p = gls(t, mag, f; σ=err)                                     # GLS power
χ2, χ0 = falc(t, mag, f; σ=err, order=2, groups=night)        # FALC χ² scan
P = 1 / f[argmin(χ2)]
fit = falc_fit(t, mag, 1 / P; σ=err, groups=night)            # amplitude, coefficients, zero points
sh = split_halves(t, mag, P; σ=err)                           # do the two halves differ?

pc = fit_phase_curve(α, reduced_magnitude(m, r, Δ); σ=err, bands=band, ref="r")
p_V, _ = geometric_albedo(pc.H, 506.0)
```

Examples (own environment `examples/Project.toml` with Plots; outputs and figures in `examples/output/<name>/`; every
example in `examples/` is also a page of the manual):
- `examples/1486_marilyn_alcdef.jl`: ALCDEF light curves from alcdef.org, light-time correction from Horizons, FALC,
  `fit_period` (P = 4.5666 ± 0.0001 h vs LCDB 4.566 h), MPB-style plot.

## Documentation, tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=docs -e 'using Pkg; Pkg.instantiate()'; julia --project=docs docs/make.jl   # → docs/build/index.html
```

## Roadmap

Observed light curves
- AoV/ANOVA periodograms; connection to AstroPDM.jl; Jurkevich method.
- Done (0.1.0-DEV): ALCDEF download and output, geometry from JPL Horizons, MPB-style plots, split-halves plots.
- ALCDEF upload helper (metadata from observation logs, cf. the old Python alctools TODO); two-body ephemeris offline.
- References: [Peranso GLS](https://www.cbabelgium.com/peranso/UserGuideHTML/Lomb-ScargleGLS.html),
  [LombScargle.jl](https://juliaastro.github.io/LombScargle.jl/stable/), [Py-PDM](https://github.com/ckm3/Py-PDM),
  [Peranso](https://www.cbabelgium.com/peranso/index.html).

Phase curves
- H–G, H–G₁G₂, H–G₁₂ phase functions and fits (H, slope); opposition effect.
- Minimum solar phase angles of TNOs for opposition-effect campaigns (`min_pha`, prototype in `src/legacy/`).

Synthetic light curves and inversion (after the v0.1.0 release, in the separate package **Russell.jl**; preparation in
work/inversion/PLAN.md, local)
- Synthetic light curves for shape models and 3D meshes; DAMIT formats.
- Light-curve inversion with the Kaasalainen method (convex inversion, native Julia, validated against the DAMIT C code).
- All-data asteroid modelling in the spirit of ADAM (light curves, occultation chords, AO images, radar).
