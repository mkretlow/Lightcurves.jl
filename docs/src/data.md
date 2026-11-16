# Data access

All functions send only the object designation (no user information). The local databases (ALCDEF archive, LCDB
release) are kept in [`datadir`](@ref) (`ENV["LIGHTCURVES_DATA"]` or the package's `data/` folder); other downloads are
cached in [`cachedir`](@ref) (`ENV["LIGHTCURVES_CACHE"]` or `~/.cache/Lightcurves.jl`).

| Source | Functions | Notes |
|---|---|---|
| MPC observations API | [`fetch_mpc_observations`](@ref), [`fetch_mpc_photometry`](@ref) | all observations in ADES form; photometry with band and error. Band corrections to V: [`mpc_to_V`](@ref) |
| JPL Horizons | [`fetch_horizons_geometry`](@ref), [`observing_geometry`](@ref) | r, Δ, phase angle, predicted V; interpolated to the observation times |
| JPL SBDB | [`fetch_sbdb`](@ref) | orbit and physical parameters (H, diameter, albedo, rotation period) |
| ALCDEF (alcdef.org) | [`alcdef_sessions`](@ref), [`fetch_alcdef`](@ref), [`read_alcdef`](@ref), [`write_alcdef`](@ref); local archive: [`fetch_alcdef_archive`](@ref), [`alcdef_archive_objects`](@ref), [`read_alcdef_archive`](@ref) | per object through the web interface, or from the full archive (≈ 135 MB, one file per object). The blocks contain the observers' contact information |
| LCDB (Warner et al.) | [`lcdb_entry`](@ref), [`fetch_lcdb_archive`](@ref), [`fetch_lcdb_member`](@ref), [`read_lcdb_summary`](@ref), [`lcdb_lookup`](@ref) | one object from its ALCDEF page, or the tables of the release zip (40 MB; local copy, or single members by HTTP range requests) |
| CdR / CourbRot files | [`read_cdr`](@ref) | R. Behrend's light-curve format |

The ALCDEF format follows the ALCDEF standard v2.3 [ALCDEF2024](@cite); the LCDB is described by [Warner2009](@citet).

```julia
using Lightcurves
ph  = fetch_mpc_photometry("612533")                       # MPC magnitudes
geo = observing_geometry("612533", ph.jd)                   # r, Δ, α from Horizons
blocks = fetch_alcdef("1486")                               # ALCDEF light curves
tab = alcdef_table(blocks)
lcdb = lcdb_entry("1486"); lcdb["Period"]
```
