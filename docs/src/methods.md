# Methods

## Period search

Times are in days (any origin; for real data use light-time corrected times), frequencies in cycles per day.

**Generalised Lomb–Scargle** ([`gls`](@ref)) fits a sinusoid plus a floating mean at each trial frequency, with weights
[ZechmeisterKurster2009](@cite). The power is the relative χ² reduction, between 0 and 1. A double-peaked asteroid
light curve gives its strongest GLS peak at *twice* the rotation frequency.

**FALC** ([`falc`](@ref), [`falc_fit`](@ref)) follows [Harris1989](@citet): at each trial frequency a Fourier series of
order ``n`` is fitted together with one zero point per group (night, station, band, …),

```math
m_i = Z_{g(i)} + \sum_{k=1}^{n} \left[a_k \cos(2\pi k f t_i) + b_k \sin(2\pi k f t_i)\right],
```

and the χ² is recorded. Order 2 or more resolves the two different halves of a double-peaked light curve, so the best
FALC frequency is the rotation frequency. The free zero points absorb calibration offsets between nights and bands.

**Spectral window** ([`spectral_window`](@ref)): ``W(f) = |\sum_i e^{2\pi i f t_i}|^2 / N^2``. Periodogram peaks at
frequencies where the window is large (one cycle per sidereal day, 1.0027 c/d, for ground-based data) are aliases.

**False-alarm probability** ([`permutation_fap`](@ref)): the data are shuffled, within groups if given, and the maximum
of the periodogram is compared with the observed one. Shuffling single points ignores correlated errors, so for sparse
survey photometry with nightly offsets the FAP is too optimistic; shuffle nightly means or pass the nights as groups.

## Folded light curves

[`fold`](@ref) gives the rotational phase, [`fourier_fit`](@ref) a weighted Fourier series in the phase.

**Split-halves test** ([`split_halves`](@ref)): the folded light curve is cut into its first and second half, which are
overlaid and binned [Harris2014](@cite). If the halves agree, the half period describes the data equally well; if they
differ, the full period is needed. A symmetric double-peaked light curve also has similar halves, so the test cannot
always decide between ``P`` and ``2P``.

## Phase curves

[`reduced_magnitude`](@ref): ``m - 5\log_{10}(r\Delta)``. [`fit_phase_curve`](@ref) fits the linear phase curve

```math
m_\mathrm{red} = H + \beta\,\alpha + \Delta_\mathrm{band},
```

with one offset per band relative to a reference band, weighted least squares and iterative σ clipping. This is
appropriate for distant objects observed at small phase angles (α ≲ 2°). [`ps1_to_V`](@ref) converts Pan-STARRS1 r and
g − r to Johnson V [Tonry2012](@cite).

## Physical quantities

[`geometric_albedo`](@ref) and [`diameter`](@ref) use ``D = 1329\,\mathrm{km}\,p_V^{-1/2}\,10^{-H_V/5}``
[HarrisLagerros2002](@cite). [`ellipsoid_amplitude`](@ref) is the light-curve amplitude of a triaxial ellipsoid at aspect
angle θ (geometric scattering, zero phase angle) [Binzel1989](@cite).

## CCD images and aperture photometry

Basic routines for reducing one's own frames to light curves; images are matrices `img[x, y]` in FITS axis order with
pixel centres at integer coordinates. [`readfits`](@ref) needs `using FITSIO` (package extension), everything else has
no extra dependency.

- **Reduction:** [`masterframe`](@ref) (pixel-wise median of bias, dark or flat frames), [`calibrate`](@ref)
  (bias, scaled dark, median-normalised flat).
- **Detection and positions:** [`findstars`](@ref) finds local maxima of the 3×3-smoothed image above a robust
  background threshold ([`robust_stats`](@ref): median and 1.4826·MAD with one 3σ clip); [`centroid`](@ref) refines
  them by iterated intensity-weighted centroids after subtracting a local ring median. It reports a status
  (`:ok`, `:edge`, `:noflux`, `:runaway` when the position moves more than `maxshift` from the start, `:maxiter`) and
  returns `NaN` coordinates on failure by default, so that a failed fit (e.g. a faint predicted position drifting onto
  a brighter neighbour) cannot be mistaken for a measurement; [`fwhm_estimate`](@ref) gives
  the median FWHM from second moments.
- **Photometry:** [`aperture`](@ref) sums a circular aperture with fractional edge pixels and subtracts the robust sky
  of an annulus; the error includes source photon noise, sky scatter and the error of the sky level
  [Howell1989](@cite). An aperture of ≈ 0.7–1 FWHM maximises the S/N for faint targets; then measure comparison stars
  with the same aperture so that the aperture correction cancels in differential photometry.
  [`instrumental_magnitude`](@ref) converts flux to magnitude.
- **Registration:** [`match_offset`](@ref) finds the translation between frames by voting on pairwise star
  differences; [`platematch`](@ref) solves an affine plate solution against catalogue tangent-plane coordinates by a
  grid search over rotation and flip at known scale, followed by a least-squares fit.

Not included: PSF fitting, image alignment by resampling, cosmic-ray rejection, WCS handling (use the header WCS or
[`platematch`](@ref)).
