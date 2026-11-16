using Lightcurves
using Test, Random, Statistics
using RecipesBase
using FITSIO                                           # loads the extension LightcurvesFITSIOExt
RecipesBase.is_key_supported(::Symbol) = true          # normally provided by Plots.jl

@testset "Lightcurves" begin

@testset "chi2_cdf" begin
    # reference values (scipy.stats.chi2.cdf)
    @test chi2_cdf(1.0, 1) ≈ 0.6826894921370859 atol = 1e-12
    @test chi2_cdf(3.84145882069412, 1) ≈ 0.95 atol = 1e-10
    for x in (0.3, 2.0, 7.5, 40.0)                       # closed forms for even k
        @test chi2_cdf(x, 2) ≈ 1 - exp(-x / 2) atol = 1e-13
        @test chi2_cdf(x, 4) ≈ 1 - exp(-x / 2) * (1 + x / 2) atol = 1e-13
    end
    @test chi2_cdf(30.0, 7) ≈ 0.99990504 atol = 2e-8         # odd-k closed form, checked by hand
    @test chi2_cdf(0.0, 3) == 0.0
    @test_throws ArgumentError chi2_cdf(1.0, 0)
end

# synthetic sparse light curve: double-peaked, P = 9.0 h, 0.2 mag, 3 points per night on 60 nights over 2 yr
rng = MersenneTwister(1)
P = 9.0 / 24
nights = sort(rand(rng, 0:730, 60))
t = sort(vcat([n .+ 0.1 .+ 0.02 .* (0:2) .+ 0.3rand(rng) for n in nights]...))
σ = fill(0.03, length(t))
lc(tt) = 0.1cos(4π * tt / P) + 0.03sin(2π * tt / P)
AMP = (m = lc.(range(0, P; length=10001)); maximum(m) - minimum(m))       # true peak-to-peak (≈ 0.235)
y = lc.(t) .+ 0.03 .* randn(rng, length(t))

@testset "frequency_grid" begin
    f = frequency_grid(t; fmin=0.5, fmax=10.0, oversample=10)
    @test f[1] == 0.5 && f[end] <= 10.0
    @test f[2] - f[1] ≈ 1 / (10 * (maximum(t) - minimum(t)))
    @test_throws ArgumentError frequency_grid([1.0, 1.0]; fmin=0.1, fmax=1.0)
end

@testset "gls" begin
    tt = collect(0:0.013:30.0); yy = sin.(2π .* tt ./ 1.7)
    f = frequency_grid(tt; fmin=0.1, fmax=2.0)
    p = gls(tt, yy, f)
    @test abs(f[argmax(p)] - 1 / 1.7) < 2e-3
    @test maximum(p) > 0.99                                   # noise-free sine (grid step limits the peak)
    @test all(0 .<= p .<= 1 + 1e-12)
    # double-peaked synthetic curve: GLS finds half the period
    f2 = frequency_grid(t; fmin=1.0, fmax=12.0, oversample=5)
    p2 = gls(t, y, f2; σ=σ)
    @test abs(f2[argmax(p2)] - 2 / P) < 0.01
end

@testset "falc" begin
    f = frequency_grid(t; fmin=1.0, fmax=12.0, oversample=5)
    χ2, χ0 = falc(t, y, f; σ=σ, order=2)
    @test χ0 > minimum(χ2)
    fb = f[argmin(χ2)]
    @test abs(fb - 1 / P) < 0.01 || abs(fb - 2 / P) < 0.01          # order 2 resolves the full period
    fit = falc_fit(t, y, 1 / P; σ=σ, order=2)
    @test fit.amplitude ≈ AMP atol = 0.02
    @test fit.chi2 / (length(t) - 5) < 2
    # zero points per group are absorbed
    g = [isodd(round(Int, ti)) ? "A" : "B" for ti in t]
    y2 = y .+ [gi == "A" ? 0.5 : -0.3 for gi in g]
    fit2 = falc_fit(t, y2, 1 / P; σ=σ, groups=g)
    @test fit2.amplitude ≈ fit.amplitude atol = 0.02
    @test sort(fit2.zero_points) ≈ sort([0.5, -0.3]) atol = 0.03
end

@testset "spectral_window and permutation_fap" begin
    w = spectral_window(t, [0.0, 1.0027, 0.37])
    @test w[1] ≈ 1.0
    @test w[2] > w[3]                                  # daily alias of nightly sampling
    f = frequency_grid(t; fmin=1.0, fmax=12.0, oversample=2)
    @test permutation_fap(gls, t, y, f; nperm=20, σ=σ, rng=MersenneTwister(2)) == 0.0
    noise = 0.03 .* randn(MersenneTwister(3), length(t))
    @test permutation_fap(gls, t, noise, f; nperm=20, rng=MersenneTwister(4)) > 0.0
end

@testset "fold, fourier_fit, split_halves" begin
    ph = fold(t, P)
    @test all(0 .<= ph .< 1)
    @test fold([0.0, 0.5, 1.0, 1.5], 1.0) == [0.0, 0.5, 0.0, 0.5]
    ff = fourier_fit(fold(t, P; t0=0.0), y; σ=σ, order=2)
    @test ff.amplitude ≈ AMP atol = 0.02
    @test ff.dchi2 > 0
    @test ff.model(0.0) ≈ ff.c0 + sum(ff.a)
    # double-peaked curve with unequal maxima: halves differ at P (small p); halves agree at 2P for a sinusoid at P
    s = split_halves(t, y, P; σ=σ, t0=0.0)
    @test s.dof > 5
    @test s.p < 0.01
    tt = collect(0:0.011:40.0); yy = 0.1sin.(2π .* tt ./ 1.3) .+ 0.01 .* randn(MersenneTwister(5), length(tt))
    s2 = split_halves(tt, yy, 2 * 1.3; σ=fill(0.01, length(tt)), t0=0.0)
    @test s2.p > 0.01
end

@testset "phase curve" begin
    @test reduced_magnitude(20.0, 10.0, 10.0) ≈ 10.0
    rng2 = MersenneTwister(6)
    α = 1.5 .* rand(rng2, 200)
    bands = rand(rng2, ["r", "g", "w"], 200)
    off = Dict("r" => 0.0, "g" => 0.55, "w" => 0.04)
    m = 4.5 .+ 0.15 .* α .+ [off[b] for b in bands] .+ 0.03 .* randn(rng2, 200)
    m[1] += 2.0                                         # an outlier
    pc = fit_phase_curve(α, m; σ=fill(0.03, 200), bands=bands, ref="r")
    @test pc.H ≈ 4.5 atol = 0.02
    @test pc.beta ≈ 0.15 atol = 0.03
    @test pc.offsets["g"] ≈ 0.55 atol = 0.02
    @test pc.offsets["r"] == 0.0
    @test !pc.keep[1] && count(pc.keep) >= 195
    pcf = fit_phase_curve(α, m; σ=fill(0.03, 200), bands=bands, ref="r", beta=0.15)
    @test pcf.beta == 0.15
    @test pcf.H ≈ 4.5 atol = 0.02
    V, eV = ps1_to_V(20.0, 0.5; g_r_err=0.05)
    @test V ≈ 20.243
    @test eV ≈ sqrt(0.0237^2 + 0.012^2)
end

@testset "physical" begin
    p, ep = geometric_albedo(4.77, 506.0; H_err=0.06, D_err=10.0)
    @test p ≈ (1329 / 506)^2 * 10^(-0.4 * 4.77)
    @test 0.08 < p < 0.09 && ep > 0
    D, _ = diameter(4.77, p)
    @test D ≈ 506.0
    @test ellipsoid_amplitude(2.0, 1.0, 1.0, π / 2) ≈ 2.5log10(2.0)
    @test ellipsoid_amplitude(2.0, 1.0, 0.8, 0.0) ≈ 0.0 atol = 1e-12
    @test 0 < ellipsoid_amplitude(2.0, 1.0, 0.8, π / 4) < 2.5log10(2.0)
    @test_throws ArgumentError ellipsoid_amplitude(1.0, 2.0, 0.5, 1.0)
end

@testset "period fit" begin
    pf = fit_period(t, y, P * 1.001; σ=σ, order=4, search=0.01)
    @test abs(pf.period - P) < 3 * pf.period_err + 1e-6
    @test pf.period_err > 0 && pf.period_err < 1e-4
    @test pf.amplitude ≈ AMP atol = 0.03
    @test 0.5 < pf.chi2_red < 2
    ph, m = phase_ephemeris(pf, [pf.tmin, pf.tmax])
    @test m[1] < m[2]                                          # minimum magnitude (brightest) before maximum
    @test std(residuals(pf)) < 0.05
    lf = lightcurve_fit(t, y, P; σ=σ, order=2)
    @test isnan(lf.period_err) && lf.period == P
    sc = fourier_order_scan(t, y, P; σ=σ, orders=1:5)
    @test argmin([s.bic for s in sc]) == 2                       # the synthetic curve has harmonics 1 and 2
    @test occursin("PeriodFit", sprint(show, pf))
end

@testset "io" begin
    @test light_time_correct(2.46e6, 1.0) ≈ 2.46e6 - 499.004783836 / 86400
    @test mpc_to_V(20.0, "R") ≈ 20.4
    @test mpc_to_V(20.0, "Pw") ≈ 19.87
    @test_throws ArgumentError mpc_to_V(20.0, "Q")
    txt = """
    STARTMETADATA
    OBJECTNUMBER=1486
    OBJECTNAME=Marilyn
    SESSIONDATE=2013-08-09
    SESSIONTIME=00:10:00
    CONTACTNAME=A. Observer
    FILTER=C
    MAGBAND=R
    EQUINOX=J2000.0
    EQUINOX=J2000.0
    DELIMITER=PIPE
    ENDMETADATA
    DATA=2456513.428390| +1.352|0.006
    DATA=2456513.429520| +1.347|
    DATA=2456513.430310| +1.367|0.006
    ENDDATA

    STARTMETADATA
    OBJECTNUMBER=1486
    SESSIONDATE=2013-08-12
    DELIMITER=SEMIC
    ENDMETADATA
    DATA=2456516.5;1.20;0.01
    DATA=2456516.6;1.25;0.01
    ENDDATA
    """
    bl = read_alcdef(IOBuffer(txt))
    @test length(bl) == 2
    @test bl[1]["OBJECTNAME"] == "Marilyn" && bl[1]["NOPE"] === nothing && bl[1]["NOPE", "x"] == "x"
    @test count(p -> p.first == "EQUINOX", bl[1].meta) == 2
    @test bl[1].jd[1] ≈ 2456513.42839 && isnan(bl[1].err[2]) && bl[2].mag[2] == 1.25
    io = IOBuffer(); write_alcdef(io, bl)
    bl2 = read_alcdef(IOBuffer(String(take!(io))))
    @test length(bl2) == 2 && bl2[2].mag ≈ bl[2].mag && bl2[2]["DELIMITER"] == "PIPE"
    noend = read_alcdef(IOBuffer("STARTMETADATA\r\nOBJECTNUMBER=1\r\nENDMETADATA\r\nDATA=1.0|2.0|0.1\r\nDATA=1.1|2.1|0.1\r\nSTARTMETADATA\r\nOBJECTNUMBER=1\r\nENDMETADATA\r\nDATA=2.0|3.0\r\n"))
    @test length(noend) == 2 && length(noend[1]) == 2 && length(noend[2]) == 1      # archive files have no ENDDATA
    tab = alcdef_table(bl)
    @test length(tab.jd) == 5 && tab.session == [1, 1, 1, 2, 2] && tab.labels[1] == "2013-08-09 A. Observer"
    f = tempname()
    write(f, "FMT DV\nNOM (1486) Marilyn\nMES A. Observer\nPOS -1 60\nTEL ROTAT\n2456520.44534   14.137 R   29.36\n2456520.44904   14.162 R   27.76\n")
    c = read_cdr(f)
    @test c.object == "(1486) Marilyn" && c.mag == [14.137, 14.162]
    @test c.jd[1] ≈ 2456520.44534 + 30 / 86400                   # start of exposure → mid-exposure
    @test c.dmag == c.mag && size(c.comps, 2) == 0
    f2 = tempname()
    write(f2, "FMT DVRR\nNOM (1069) Planckia\nPOS 0 0\n2455684.35423  -9.29  -11.548  -9.195\n#2455684.36   -9.0  -11.0  -9.0\n2455684.36000  -9.30  -11.500  -9.200\n")
    c2 = read_cdr(f2)
    @test length(c2.jd) == 2 && size(c2.comps) == (2, 2)                  # the commented line is skipped
    @test c2.dmag[1] ≈ -9.29 - (-11.548 - 9.195) / 2
end

@testset "parsers of fetched data (offline)" begin
    htxt = "header\n\$\$SOE\n 2026-Nov-01 00:00, , ,  20.565,   n.a.,   37.64549, -0.5,  36.85460, -20.1,    1.4253,\n 2026-Nov-02 00:00, , ,  20.563,   n.a.,   37.64540, -0.5,  36.84300, -20.0,    1.4160,\n\$\$EOE\nfooter"
    g = parse_horizons_geometry(htxt)
    @test length(g.jd) == 2 && g.jd[1] ≈ 2461345.5 && g.r[1] ≈ 37.64549 && g.alpha[2] ≈ 1.416 && g.vmag[1] ≈ 20.565
    @test Lightcurves._iso_to_jd("2026-11-11T05:55:00.000Z") ≈ 2461355.5 + (5 + 55 / 60) / 24
    f = tempname()
    write(f, """
    ASTEROID LIGHTCURVE DATABASE (LCDB) SUMMARY (Complete)

    NUMBER    NAME                           DESIG                FAM      S CLASS      S F DIA.     S H      B  S G      G1     G2     S F ALBEDO F PERIOD        P DESC          F AMIN AMAX U  NOTES BIN SAM SurvA NEX PRI
    -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     136108   Haumea                         2003 EL61            9301     A C*         L   1237.24  L  0.425    A  0.12                L   0.78        3.9154                       0.05 0.32 3         B
       1486   Marilyn                        1938 QA              9104     A S*         C      7.96  M 12.86     A  0.24                A   0.20        4.566                        0.42 0.48 3
    """)
    tab = read_lcdb_summary(f)
    @test tab.number == [136108, 1486] && tab.name == ["Haumea", "Marilyn"]
    @test tab.period ≈ [3.9154, 4.566] && tab.albedo ≈ [0.78, 0.20] && tab.U == ["3", "3"] && tab.amax ≈ [0.32, 0.48]
    @test lcdb_lookup(tab, 1486) == 2 && lcdb_lookup(tab, "2003 EL61") == 1 && lcdb_lookup(tab, "nope") === nothing
    @test Lightcurves._urlencode("a b;'é") == "a+b%3B%27%C3%A9"
end

@testset "photometry (images)" begin
    # synthetic frame: Gaussian stars on a flat background with Gaussian noise
    g(x, y, x0, y0, F, σ) = F / (2π * σ^2) * exp(-((x - x0)^2 + (y - y0)^2) / (2σ^2))
    rng = MersenneTwister(7)
    σpsf = 1.5; fw = 2.3548 * σpsf
    stars = [(40.3, 50.7, 2e5), (120.6, 80.2, 1e5), (70.1, 140.9, 5e4), (150.4, 150.5, 3e4), (30.8, 160.2, 2e4)]
    model(dx=0.0, dy=0.0) = [100.0 + sum(g(i, j, x + dx, y + dy, F, σpsf) for (x, y, F) in stars) for i in 1:200, j in 1:200]
    noise = 5.0
    img = model() .+ noise .* randn(rng, 200, 200)

    lv, σbg = robust_stats(img)
    @test isapprox(lv, 100; atol=0.3) && isapprox(σbg, noise; rtol=0.05)
    @test robust_stats(fill(3.0, 10)) == (3.0, 0.0)

    st = findstars(img; nsigma=5, fwhm=3)
    @test length(st) == length(stars)
    @test all(s -> any(p -> hypot(s.x - p[1], s.y - p[2]) <= 1, stars), st)
    @test issorted([s.peak for s in st]; rev=true)
    c = centroid(img, 40, 51)
    @test hypot(c.x - 40.3, c.y - 50.7) < 0.02 && c.status === :ok
    xc, yc = centroid(img, 40, 51); @test (xc, yc) == (c.x, c.y)   # destructuring ignores the status
    # failures: NaN position by default, start or last iterate on request, status always reported
    e = centroid(img, 2, 2)                                        # box outside the image
    @test e.status === :edge && isnan(e.x) && isnan(e.y)
    @test centroid(img, 2, 2; onfail=:start) == (x=2.0, y=2.0, status=:edge)
    @test centroid(img, NaN, 50).status === :edge
    @test centroid(fill(100.0, 50, 50), 25, 25).status === :noflux # blank, noiseless sky
    w = centroid(img, 45.5, 50.7)                                  # 5.2 px from a bright star: drifts onto it
    @test w.status === :runaway && isnan(w.x)
    wl = centroid(img, 45.5, 50.7; onfail=:last)
    @test wl.status === :runaway && hypot(wl.x - 45.5, wl.y - 50.7) > 4
    wm = centroid(img, 45.5, 50.7; maxshift=7)                     # allowed to move that far: finds the star
    @test wm.status === :ok && hypot(wm.x - 40.3, wm.y - 50.7) < 0.05
    @test centroid(img, 40, 51; iter=1).status === :maxiter        # one step from 0.7 px away does not converge
    @test_throws ArgumentError centroid(img, 40, 51; onfail=:zero)
    sat = min.(img, 100.0 + 0.3 * maximum(img .- 100.0))           # flat-topped (saturated) star
    cs = centroid(sat, 40, 51)
    @test cs.status === :ok && hypot(cs.x - 40.3, cs.y - 50.7) < 0.1

    a = aperture(img, c.x, c.y; r=4fw, rin=5fw, rout=7fw)
    @test isapprox(a.flux, 2e5; rtol=0.01)
    @test abs(a.flux - 2e5) < 3a.σflux
    @test isapprox(a.npix, π * (4fw)^2; rtol=0.01)
    @test isapprox(a.sky, 100; atol=0.5)
    @test isnan(aperture(img, 3, 100; r=4).flux)                  # aperture leaves the image
    @test isnan(aperture(img, NaN, NaN; r=4).flux)                # failed centroid passed on
    m = instrumental_magnitude(a.flux, a.σflux; zp=25)
    @test m.mag ≈ 25 - 2.5log10(a.flux)
    @test isapprox(m.σmag, 1.0857362 * a.σflux / a.flux; rtol=1e-6)
    @test isnan(instrumental_magnitude(-1.0).mag)

    @test isapprox(fwhm_estimate(img, st; r=6), fw; rtol=0.1)

    img2 = model(7.4, -3.2) .+ noise .* randn(rng, 200, 200)
    p1 = [centroid(img, s.x, s.y) for s in st]
    p2 = [centroid(img2, s.x, s.y) for s in findstars(img2)]
    o = match_offset(p1, p2; tol=1.0)
    @test o.nmatch == length(stars) && isapprox(o.dx, 7.4; atol=0.1) && isapprox(o.dy, -3.2; atol=0.1)

    # plate match: sky = A·pix + t with scale 0.5"/px, rotation 30°, flipped; extra catalogue stars
    A = 0.5 .* [cosd(30) -sind(30); sind(30) cosd(30)] * [-1 0; 0 1]; t0 = [12.0, -40.0]
    pix = [(p.x, p.y) for p in p1]
    sky = [Tuple(A * [p[1], p[2]] .+ t0) for p in pix]
    sky = vcat(sky[1:2], [(500.0, 500.0), (-300.0, 20.0)], sky[3:end])
    pm = platematch(pix, sky; scale=0.5, tol=0.5)
    @test pm.nmatch == length(pix) && isapprox(pm.A, A; atol=1e-3) && isapprox(pm.t, t0; atol=0.1)

    # reduction
    bias = fill(300.0, 20, 20); flat = 2000 .* (1 .+ 0.1 .* (1:20) ./ 20) .* ones(20)'
    raw = bias .+ 50 .* flat ./ median(vec(flat))
    @test calibrate(raw; bias=bias, flat=flat) ≈ fill(50.0, 20, 20)
    @test masterframe([fill(1.0, 3, 3), fill(5.0, 3, 3), fill(2.0, 3, 3)]) == fill(2.0, 3, 3)
    @test_throws DimensionMismatch masterframe([ones(2, 2), ones(3, 3)])

    # FITS through the FITSIO extension
    f = joinpath(mktempdir(), "t.fits")
    let
        FITS(f, "w") do io
            write(io, Float32.(img[1:50, 1:60]); header=FITSHeader(["DATE-OBS"], ["2025-09-23T01:02:03"], [""]))
        end
        im, h = readfits(f)
        @test size(im) == (50, 60) && isapprox(im, img[1:50, 1:60]; rtol=1e-6) && h["DATE-OBS"] == "2025-09-23T01:02:03"
        @test masterframe([f, f]) ≈ im
    end
end

@testset "recipes" begin
    pf = lightcurve_fit(t, y, P; σ=σ, order=2, groups=[ti < 365 ? "Y1" : "Y2" for ti in t], labels=["2024", "2025"])
    series = RecipesBase.apply_recipe(Dict{Symbol,Any}(), pf)
    @test length(series) >= 3                                    # two groups + Fourier curve
    @test any(s -> get(s.plotattributes, :label, "") == "Fourier order 2", series)
    bl = read_alcdef(IOBuffer("STARTMETADATA\nSESSIONDATE=2020-01-01\nENDMETADATA\nDATA=2458850.5|1.0|0.01\nDATA=2458850.6|1.1|0.01\nENDDATA\n"))
    @test length(RecipesBase.apply_recipe(Dict{Symbol,Any}(), bl)) == 1
end

if get(ENV, "LIGHTCURVES_NETWORK_TESTS", "false") == "true"
    @testset "network (opt-in)" begin
        ph = fetch_mpc_photometry("612533"); @test length(ph.jd) > 100
        @test fetch_sbdb("612533")["phys"]["H"].value !== nothing
        @test length(alcdef_sessions("1486")) >= 5
        @test get(lcdb_entry("1486"), "Period", "") == "4.566"
    end
end

end
