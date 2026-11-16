# CCD reduction, star detection, centroids, aperture photometry and simple plate matching.
#
# Taken over from photlib.jl of the 2010 SS43 occultation workspace (OSN 1.5 m and LCO 0.4 m frames; user, 2026-11-16).
# Images are matrices img[x, y] in FITS axis order; pixel centres lie at integer coordinates (1-based), so the
# centre of the first pixel is (1, 1). Reading FITS files needs FITSIO.jl (package extension LightcurvesFITSIOExt).

# set by the FITSIO extension
const _FITS_READER = Ref{Any}(nothing)

"""
    readfits(file; hdu=nothing) -> (img::Matrix{Float64}, header)

Read a two-dimensional FITS image as `img[x, y]` (FITS axis order) with its header. `hdu` selects the HDU (1-based);
by default the first HDU that holds a 2-D image, so tile-compressed (fpack) files work too. The header is a
`FITSIO.FITSHeader` (`header["DATE-OBS"]`).

Needs FITSIO.jl: `using FITSIO` loads the package extension.
"""
function readfits(file::AbstractString; hdu=nothing)
    _FITS_READER[] === nothing && error("readfits needs FITSIO.jl: run `using FITSIO` (package extension of Lightcurves)")
    return _FITS_READER[](file, hdu)
end

"""
    robust_stats(v) -> (level, σ)

Robust level and scatter of the values `v`: median and 1.4826·MAD, recomputed once after clipping at 3σ. Used for the
image background and the sky annulus.
"""
function robust_stats(v)
    v = filter(isfinite, vec(v))
    isempty(v) && return (NaN, NaN)
    m = median(v); s = 1.4826 * median(abs.(v .- m))
    s > 0 || return (m, 0.0)
    k = abs.(v .- m) .< 3s
    m = median(v[k]); s = 1.4826 * median(abs.(v[k] .- m))
    return (m, s)
end

"""
    masterframe(frames) -> Matrix{Float64}

Pixel-wise median of a stack of frames (bias, dark or flat). `frames` is a vector of equally sized matrices or of FITS
file names (the latter needs FITSIO.jl, see [`readfits`](@ref)).
"""
function masterframe(frames::AbstractVector{<:AbstractMatrix})
    isempty(frames) && throw(ArgumentError("no frames"))
    sz = size(first(frames))
    all(f -> size(f) == sz, frames) || throw(DimensionMismatch("frames differ in size"))
    out = Matrix{Float64}(undef, sz)
    buf = Vector{Float64}(undef, length(frames))
    for I in CartesianIndices(out)
        for (k, f) in enumerate(frames)
            buf[k] = f[I]
        end
        out[I] = median!(buf)
    end
    return out
end
masterframe(files::AbstractVector{<:AbstractString}) = masterframe([readfits(f)[1] for f in files])

"""
    calibrate(img; bias=nothing, dark=nothing, flat=nothing) -> Matrix{Float64}

Basic CCD reduction ``(\\mathrm{img} - \\mathrm{bias} - \\mathrm{dark}) / (\\mathrm{flat}/\\mathrm{median(flat)})``.
`dark` is a bias-subtracted dark scaled to the exposure time of `img`; `flat` a bias- (and dark-) subtracted flat
field, normalised here by its median. Each frame may be omitted.
"""
function calibrate(img::AbstractMatrix; bias=nothing, dark=nothing, flat=nothing)
    out = Float64.(img)
    bias === nothing || (out .-= bias)
    dark === nothing || (out .-= dark)
    if flat !== nothing
        f = flat ./ median(vec(flat))
        out ./= f
    end
    return out
end

"""
    findstars(img; nsigma=5.0, fwhm=3.0, border=8) -> Vector{(x, y, peak)}

Star detection: local maxima of the 3×3 box-smoothed image that lie `nsigma`·σ (σ of the 3×3 mean) above the robust
background ([`robust_stats`](@ref)) and are separated by at least `fwhm` pixels; pixels within `border` of the edge are
skipped. Returns named tuples `(x, y, peak)` (integer pixel positions, peak above background in the smoothed image),
brightest first. Refine the positions with [`centroid`](@ref).
"""
function findstars(img::AbstractMatrix; nsigma=5.0, fwhm=3.0, border=8)
    nx, ny = size(img)
    bg, σ = robust_stats(img)
    sm = fill(bg, nx, ny)
    for j in 2:ny-1, i in 2:nx-1
        sm[i, j] = sum(@view img[i-1:i+1, j-1:j+1]) / 9
    end
    thr = bg + nsigma * σ / 3
    r = max(1, round(Int, fwhm))
    b = max(border, r + 1)
    cand = NamedTuple{(:x, :y, :peak),NTuple{3,Float64}}[]
    for j in b+1:ny-b, i in b+1:nx-b
        v = sm[i, j]
        v > thr || continue
        v == maximum(@view sm[i-r:i+r, j-r:j+r]) || continue
        push!(cand, (x=i, y=j, peak=v - bg))
    end
    sort!(cand; by=s -> -s.peak)
    # flat-topped (saturated) stars give several equal maxima: keep the first within r
    out = similar(cand, 0)
    for s in cand
        any(o -> hypot(o.x - s.x, o.y - s.y) <= r, out) || push!(out, s)
    end
    return out
end

"""
    centroid(img, x, y; r=4, iter=10, tol=0.01, maxshift=r, onfail=:nan) -> (x, y, status)

Intensity-weighted centroid within radius `r` around the current position, after subtracting the median of a ring
(r + 2 … r + 6 px), iterated until the shift per step is below `tol` px in both coordinates.

`status` tells whether the fit succeeded:

- `:ok`: converged within `maxshift` of the start position;
- `:edge`: the box (radius r + 6) leaves the image, or the start position is not finite;
- `:noflux`: no pixel above the ring median (blank sky);
- `:runaway`: the centroid moved more than `maxshift` px from the start position (in total, not per step), e.g. a
  faint or empty position drifting onto a brighter neighbour;
- `:maxiter`: not converged after `iter` iterations.

On failure the returned position is `NaN` by default (`onfail=:nan`), so that an unchecked failure cannot pass as a
measured position; `onfail=:start` returns the start position and `onfail=:last` the last iterate instead (the
`status` still reports the failure). `xc, yc = centroid(...)` ignores the status.
"""
function centroid(img::AbstractMatrix, x, y; r=4, iter=10, tol=0.01, maxshift=r, onfail::Symbol=:nan)
    onfail in (:nan, :start, :last) || throw(ArgumentError("onfail must be :nan, :start or :last, got :$onfail"))
    nx, ny = size(img)
    x, y = Float64(x), Float64(y)
    x0, y0 = x, y
    R = ceil(Int, r)
    function fail(status, xl, yl)
        onfail === :nan && return (x=NaN, y=NaN, status=status)
        onfail === :start && return (x=x0, y=y0, status=status)
        return (x=xl, y=yl, status=status)
    end
    (isfinite(x) && isfinite(y)) || return fail(:edge, x, y)
    for _ in 1:iter
        i0, j0 = round(Int, x), round(Int, y)
        (i0 - R - 6 < 1 || j0 - R - 6 < 1 || i0 + R + 6 > nx || j0 + R + 6 > ny) && return fail(:edge, x, y)
        ring = Float64[]
        for j in j0-R-6:j0+R+6, i in i0-R-6:i0+R+6
            d = hypot(i - x, j - y)
            r + 2 < d < r + 6 && push!(ring, img[i, j])
        end
        b = median(ring)
        sw = sx = sy = 0.0
        for j in j0-R:j0+R, i in i0-R:i0+R
            hypot(i - x, j - y) <= r || continue
            w = max(img[i, j] - b, 0.0)
            sw += w; sx += w * i; sy += w * j
        end
        sw > 0 || return fail(:noflux, x, y)
        xn, yn = sx / sw, sy / sw
        hypot(xn - x0, yn - y0) > maxshift && return fail(:runaway, xn, yn)
        (abs(xn - x) < tol && abs(yn - y) < tol) && return (x=xn, y=yn, status=:ok)
        x, y = xn, yn
    end
    return fail(:maxiter, x, y)
end

"""
    aperture(img, x, y; r=4.0, rin=8.0, rout=13.0, gain=1.0)
        -> (flux, σflux, sky, σsky, npix)

Circular aperture photometry at `(x, y)` with radius `r` and a sky annulus `rin … rout` (pixels). Pixels cut by the
aperture edge are weighted by their covered fraction (5×5 sub-pixel sampling). The sky level is the robust median of
the annulus ([`robust_stats`](@ref)). Flux and sky in ADU; `npix` the effective aperture area.

The error combines the photon noise of the source (`gain` in e⁻/ADU), the measured sky scatter per pixel (which
includes read noise) and the error of the mean sky level [Howell1989](@cite):
``σ^2 = F/g + n_\\mathrm{pix}σ_\\mathrm{sky}^2 + n_\\mathrm{pix}^2σ_\\mathrm{sky}^2/n_\\mathrm{sky}``.
Returns `NaN` flux if the position is not finite (e.g. a failed [`centroid`](@ref)), the aperture leaves the image
or the annulus has fewer than 10 pixels.
"""
function aperture(img::AbstractMatrix, x, y; r=4.0, rin=8.0, rout=13.0, gain=1.0)
    nx, ny = size(img)
    nan = (flux=NaN, σflux=NaN, sky=NaN, σsky=NaN, npix=NaN)
    (isfinite(x) && isfinite(y)) || return nan                  # e.g. a failed centroid
    (x - r < 0.5 || y - r < 0.5 || x + r > nx + 0.5 || y + r > ny + 0.5) && return nan
    R = ceil(Int, rout) + 1
    i0, j0 = round(Int, x), round(Int, y)
    sky = Float64[]; f = 0.0; npix = 0.0
    for j in max(1, j0 - R):min(ny, j0 + R), i in max(1, i0 - R):min(nx, i0 + R)
        d = hypot(i - x, j - y)
        rin <= d <= rout && push!(sky, img[i, j])
        if d <= r - 0.75
            f += img[i, j]; npix += 1
        elseif d <= r + 0.75
            w = 0.0
            for a in -0.4:0.2:0.4, b in -0.4:0.2:0.4
                w += hypot(i + a - x, j + b - y) <= r
            end
            w /= 25
            f += w * img[i, j]; npix += w
        end
    end
    length(sky) < 10 && return nan
    s, σs = robust_stats(sky)
    flux = f - npix * s
    var = max(flux, 0.0) / gain + npix * σs^2 + npix^2 * σs^2 / length(sky)
    return (flux=flux, σflux=sqrt(var), sky=s, σsky=σs, npix=npix)
end

"""
    instrumental_magnitude(flux, σflux=0; zp=0.0) -> (mag, σmag)

Magnitude ``zp - 2.5\\log_{10} F`` and its error ``1.0857\\,σ_F/F``. Non-positive fluxes give `NaN`.
"""
function instrumental_magnitude(flux::Real, σflux::Real=0; zp=0.0)
    flux > 0 || return (mag=NaN, σmag=NaN)
    return (mag=zp - 2.5 * log10(flux), σmag=2.5 / log(10) * σflux / flux)
end

"""
    fwhm_estimate(img, stars; r=6, n=15) -> Float64

Median FWHM [px] of the `n` brightest `stars` (output of [`findstars`](@ref) or `(x, y)` tuples) from the second
moments within radius `r` after background subtraction, assuming circular Gaussian profiles
(FWHM = 2.3548 σ). The moments are truncated at `r`, so choose `r` ≳ 1.5 FWHM.
"""
function fwhm_estimate(img::AbstractMatrix, stars; r=6, n=15)
    w = Float64[]
    for s in stars[1:min(n, length(stars))]
        xc, yc, st = centroid(img, s[1], s[2]; r=r)
        st === :ok || continue
        i0, j0 = round(Int, xc), round(Int, yc)
        (i0 - r - 8 < 1 || j0 - r - 8 < 1 || i0 + r + 8 > size(img, 1) || j0 + r + 8 > size(img, 2)) && continue
        ring = [img[i, j] for j in j0-r-8:j0+r+8, i in i0-r-8:i0+r+8 if r + 3 < hypot(i - xc, j - yc) < r + 8]
        b = median(ring)
        s0 = sxx = syy = 0.0
        for j in j0-r:j0+r, i in i0-r:i0+r
            v = img[i, j] - b
            v > 0 || continue
            s0 += v; sxx += v * (i - xc)^2; syy += v * (j - yc)^2
        end
        s0 > 0 && push!(w, 2.3548 * sqrt((sxx + syy) / (2s0)))
    end
    return isempty(w) ? NaN : median(w)
end

"""
    match_offset(a, b; tol=2.0, nmax=25) -> (dx, dy, nmatch)

Translation `b − a` between two star lists (`(x, y)` tuples, brightest first; only the first `nmax` are used): every
pairwise difference is tried and the one under which most stars of `a` have a partner in `b` within `tol` pixels wins.
For frames of one night without rotation; use [`platematch`](@ref) for rotation and scale.
"""
function match_offset(a, b; tol=2.0, nmax=25)
    a = a[1:min(nmax, length(a))]; b = b[1:min(nmax, length(b))]
    best = (dx=0.0, dy=0.0, nmatch=0)
    for p in a, q in b
        dx, dy = q[1] - p[1], q[2] - p[2]
        n = count(pp -> any(qq -> hypot(qq[1] - pp[1] - dx, qq[2] - pp[2] - dy) < tol, b), a)
        n > best.nmatch && (best = (dx=float(dx), dy=float(dy), nmatch=n))
    end
    return best
end

"""
    platematch(pix, sky; scale, angles=0:1:359, flips=(false, true), tol=2.0, npix=20, nsky=60)
        -> (A, t, nmatch, pairs)

Affine plate solution ``(ξ, η) ≈ A\\,(x, y) + t`` from image stars `pix` [(x, y), …] and catalogue tangent-plane
offsets `sky` [(ξ, η) in arcsec, …], both brightest first. Brute-force search over rotation angle (degrees) and flip at
the given `scale` [arcsec/px] with translation voting (`tol` in arcsec), then a linear least-squares fit of the six
affine parameters to all matched pairs (if ≥ 3). `pairs` lists the matched indices `(i_pix, j_sky)`.
"""
function platematch(pix, sky; scale, angles=0:1:359, flips=(false, true), tol=2.0, npix=20, nsky=60)
    P = pix[1:min(npix, length(pix))]; Q = sky[1:min(nsky, length(sky))]
    best = (n=0, R=zeros(2, 2), t=zeros(2))
    for fl in flips, θ in angles
        R = scale .* [cosd(θ) -sind(θ); sind(θ) cosd(θ)] * (fl ? [-1 0; 0 1] : [1 0; 0 1])
        Pr = [R * [p[1], p[2]] for p in P]
        for p in Pr[1:min(8, end)], q in Q
            t = [q[1], q[2]] .- p
            n = count(pp -> any(qq -> hypot(qq[1] - pp[1] - t[1], qq[2] - pp[2] - t[2]) < tol, Q), Pr)
            n > best.n && (best = (n=n, R=R, t=t))
        end
    end
    pairs = Tuple{Int,Int}[]
    for (i, p) in enumerate(pix)
        s = best.R * [p[1], p[2]] .+ best.t
        d, j = findmin([hypot(q[1] - s[1], q[2] - s[2]) for q in sky])
        d < tol && push!(pairs, (i, j))
    end
    if length(pairs) >= 3
        X = permutedims(hcat([[pix[i][1], pix[i][2], 1.0] for (i, _) in pairs]...))
        Y = permutedims(hcat([[sky[j][1], sky[j][2]] for (_, j) in pairs]...))
        C = X \ Y
        return (A=Matrix(permutedims(C[1:2, :])), t=C[3, :], nmatch=length(pairs), pairs=pairs)
    end
    return (A=best.R, t=best.t, nmatch=best.n, pairs=pairs)
end
