# API

```@docs
Lightcurves
```

```@index
```

## Period search

```@docs
frequency_grid
gls
falc
falc_fit
spectral_window
permutation_fap
```

## Folded light curves

```@docs
fold
fourier_fit
split_halves
```

## Phase curves

```@docs
reduced_magnitude
fit_phase_curve
ps1_to_V
```

## Physical quantities

```@docs
geometric_albedo
diameter
ellipsoid_amplitude
```

## Period fit and light-curve model

```@docs
PeriodFit
fit_period
lightcurve_fit
fourier_series
fourier_order_scan
phase_ephemeris
residuals
```

## File formats and corrections

```@docs
ALCDEFBlock
read_alcdef
write_alcdef
alcdef_table
read_cdr
light_time_correct
MPC_BAND_TO_V
mpc_to_V
```

## Data access

```@docs
fetch_mpc_observations
fetch_mpc_photometry
fetch_horizons_geometry
parse_horizons_geometry
observing_geometry
fetch_sbdb
fetch_alcdef
alcdef_sessions
lcdb_entry
fetch_lcdb_member
read_lcdb_summary
lcdb_lookup
cachedir
datadir
fetch_lcdb_archive
fetch_alcdef_archive
alcdef_archive_objects
read_alcdef_archive
```

## CCD images and aperture photometry

```@docs
readfits
robust_stats
masterframe
calibrate
findstars
centroid
aperture
instrumental_magnitude
fwhm_estimate
match_offset
platematch
```

## Plot recipes (Plots.jl)

`plot(fit::PeriodFit; cycles=1, errorbars=true, residuals=false)` draws the phased light curve in the style of the Minor
Planet Bulletin; `plot(blocks::Vector{ALCDEFBlock}; reltime=true)` the raw ALCDEF sessions.

```@docs
periodogramplot
periodogramplot!
phasecurveplot
phasecurveplot!
splithalvesplot
splithalvesplot!
windowplot
windowplot!
```

## Statistics

```@docs
chi2_cdf
```
