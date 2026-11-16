"""
    Lightcurves

Rotational light curves and phase curves of small Solar System bodies (asteroids, Centaurs, TNOs).

- Period search: generalised Lomb–Scargle ([`gls`](@ref)), FALC Fourier-series χ² scan ([`falc`](@ref), [`falc_fit`](@ref)),
  spectral window ([`spectral_window`](@ref)), frequency grids ([`frequency_grid`](@ref)), permutation false-alarm
  probabilities ([`permutation_fap`](@ref)).
- Folding: [`fold`](@ref), Fourier fits ([`fourier_fit`](@ref)), split-halves test ([`split_halves`](@ref)).
- Phase curves: reduced magnitudes ([`reduced_magnitude`](@ref)), linear phase-curve fit with band offsets and clipping
  ([`fit_phase_curve`](@ref)), photometric transformations ([`ps1_to_V`](@ref)).
- Physical quantities: geometric albedo and diameter ([`geometric_albedo`](@ref), [`diameter`](@ref)), light-curve
  amplitude of a triaxial ellipsoid ([`ellipsoid_amplitude`](@ref)).

- Period fit with the period free ([`fit_period`](@ref), result [`PeriodFit`](@ref)), Fourier-order scan, phase ephemeris.
- File formats: ALCDEF ([`read_alcdef`](@ref), [`write_alcdef`](@ref)), CdR / CourbRot ([`read_cdr`](@ref)); MPC band
  corrections ([`mpc_to_V`](@ref)); light-time correction.
- Data access: MPC observations ([`fetch_mpc_photometry`](@ref)), JPL Horizons geometry ([`observing_geometry`](@ref)),
  JPL SBDB ([`fetch_sbdb`](@ref)), ALCDEF light curves ([`fetch_alcdef`](@ref)), LCDB ([`lcdb_entry`](@ref),
  [`fetch_lcdb_member`](@ref), [`read_lcdb_summary`](@ref)).
- CCD images: reduction ([`masterframe`](@ref), [`calibrate`](@ref)), star detection ([`findstars`](@ref)), centroids,
  aperture photometry ([`aperture`](@ref), [`instrumental_magnitude`](@ref)), FWHM, frame offsets and a simple plate
  match ([`platematch`](@ref)); FITS reading ([`readfits`](@ref)) with `using FITSIO` (package extension).
- Plot recipes for Plots.jl: `plot(::PeriodFit)` (phased light curve, Minor Planet Bulletin style),
  `plot(::Vector{ALCDEFBlock})`, [`periodogramplot`](@ref), [`phasecurveplot`](@ref), [`splithalvesplot`](@ref),
  [`windowplot`](@ref).

Light dependencies only (RecipesBase, JSON, CodecZlib and standard libraries). PDM: package AstroPDM.jl.

Author: Mike Kretlow <mike@kretlow.de>
"""
module Lightcurves

using LinearAlgebra
using Printf
using Random
using Statistics
using Dates: Dates, DateTime
using Downloads
using JSON
using CodecZlib
using RecipesBase

include("stats.jl")
include("periodograms.jl")
include("folding.jl")
include("phasecurves.jl")
include("physical.jl")
include("io.jl")
include("periodfit.jl")
include("fetch.jl")
include("photometry.jl")
include("recipes.jl")

export gls, falc, falc_fit, spectral_window, frequency_grid, permutation_fap
export fold, fourier_fit, split_halves
export reduced_magnitude, fit_phase_curve, ps1_to_V
export geometric_albedo, diameter, ellipsoid_amplitude
export chi2_cdf
export PeriodFit, fit_period, lightcurve_fit, fourier_series, fourier_order_scan, phase_ephemeris, residuals
export ALCDEFBlock, read_alcdef, write_alcdef, alcdef_table, read_cdr, light_time_correct, MPC_BAND_TO_V, mpc_to_V
export fetch_mpc_observations, fetch_mpc_photometry, fetch_horizons_geometry, parse_horizons_geometry,
       observing_geometry, fetch_sbdb, fetch_alcdef, alcdef_sessions, lcdb_entry, fetch_lcdb_member,
       read_lcdb_summary, lcdb_lookup, cachedir, datadir, fetch_lcdb_archive, fetch_alcdef_archive, alcdef_archive_objects,
       read_alcdef_archive
export readfits, robust_stats, masterframe, calibrate, findstars, centroid, aperture, instrumental_magnitude,
       fwhm_estimate, match_offset, platematch
export periodogramplot, periodogramplot!, phasecurveplot, phasecurveplot!, splithalvesplot, splithalvesplot!,
       windowplot, windowplot!

end # module
