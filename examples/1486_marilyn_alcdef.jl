# (1486) Marilyn: rotation period from ALCDEF light curves, fetched from alcdef.org, compared with the LCDB.
#
#   julia --project=examples examples/1486_marilyn_alcdef.jl > examples/output/1486_marilyn_alcdef/output.txt
#
# Needs network access the first time: the ALCDEF blocks are cached in cachedir() (they contain the observers' contact
# information and are therefore not stored in the repository), the geometry comes from JPL Horizons.
# Figures in examples/output/1486_marilyn_alcdef/.

using Lightcurves, Printf, Statistics, Plots
ENV["GKSwstype"] = "100"
default(fontfamily="sans-serif", dpi=150, grid=true, gridalpha=0.15)
const OUT = mkpath(joinpath(@__DIR__, "output", "1486_marilyn_alcdef"))
section(s) = println("\n", "="^78, "\n", s, "\n", "="^78)

section("1. Data: ALCDEF light curves and the LCDB entry")
cache = joinpath(cachedir(), "alcdef_1486.txt")
if !isfile(cache)
    mkpath(dirname(cache)); write_alcdef(cache, fetch_alcdef("1486"))
end
blocks = read_alcdef(cache)
foreach(println, blocks)
lcdb = lcdb_entry("1486")
@printf("LCDB: P = %s h, amplitude %s–%s mag, U = %s\n", lcdb["Period"], lcdb["AMin"], lcdb["AMax"], lcdb["U"])
savefig(plot(blocks; title="(1486) Marilyn: ALCDEF sessions", size=(1000, 450), left_margin=5Plots.mm, bottom_margin=5Plots.mm),
        joinpath(OUT, "sessions_raw.png"))

section("2. Light-time correction (geocentric distance from JPL Horizons)")
tab = alcdef_table(blocks)
geo = observing_geometry("1486", tab.jd; step="6h")
t = light_time_correct(tab.jd, geo.delta)
@printf("Δ = %.3f–%.3f au → light time %.1f–%.1f min; phase angle %.1f–%.1f°\n", extrema(geo.delta)...,
        24 * 60 * Lightcurves.AU_LIGHT_DAYS .* extrema(geo.delta)..., extrema(geo.alpha)...)

section("3. Period search (FALC with one zero point per session; GLS)")
σ = [isnan(e) ? 0.01 : e for e in tab.err]
f = frequency_grid(t; fmin=24 / 20, fmax=24 / 2, oversample=10)
χ2, χ0 = falc(t, tab.mag, f; σ=σ, order=4, groups=tab.session)
P = gls(t, tab.mag .- [mean(tab.mag[tab.session .== s]) for s in tab.session], f; σ=σ)
Pf = 24 / f[argmin(χ2)]
@printf("FALC (order 4) best P = %.4f h; GLS best P = %.4f h (half the rotation period for a double-peaked curve)\n",
        Pf, 24 / f[argmax(P)])
savefig(plot(periodogramplot(f, χ0 .- χ2; kind=:falc, periods=[Pf, Pf / 2], title="FALC, order 4"),
             periodogramplot(f, P; kind=:gls, title="GLS"); layout=(2, 1), size=(1000, 650), left_margin=5Plots.mm),
        joinpath(OUT, "periodograms.png"))

section("4. Period fit with the period free; Fourier order")
for r in fourier_order_scan(t, tab.mag, Pf / 24; σ=σ, groups=tab.session, orders=2:8)
    @printf("order %d: χ²_red %.2f, amplitude %.3f, BIC %.1f\n", r.order, r.chi2_red, r.amplitude, r.bic)
end
pf = fit_period(t, tab.mag, Pf / 24; σ=σ, order=4, groups=tab.session, labels=tab.labels, object="(1486) Marilyn")
println(pf)
@printf("P = %.4f ± %.4f h (LCDB %s h), amplitude %.2f mag, epoch of maximum brightness JD %.5f\n", 24pf.period,
        24pf.period_err, lcdb["Period"], pf.amplitude, pf.tmin)
savefig(plot(pf; residuals=true, left_margin=5Plots.mm), joinpath(OUT, "phased_lightcurve_mpb.png"))
savefig(plot(pf; cycles=2, errorbars=false, left_margin=5Plots.mm), joinpath(OUT, "phased_lightcurve_2cycles.png"))
sh = split_halves(t, tab.mag .- pf.zero_points[pf.g], pf.period; σ=σ)
@printf("split halves at P: χ² = %.1f for %d bins, p = %.3g (two different halves: double-peaked)\n", sh.chi2, sh.dof, sh.p)
savefig(splithalvesplot(t, tab.mag .- pf.zero_points[pf.g], pf.period; σ=σ, size=(900, 450)), joinpath(OUT, "split_halves.png"))
println("\nfigures: ", join(sort(readdir(OUT)), ", "))
