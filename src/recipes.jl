# Plot recipes (RecipesBase; used by Plots.jl when it is loaded). Styles follow the plots of the 2002 XV93 analysis and
# the phased light curves of the Minor Planet Bulletin (Δmag inverted, one colour per session, Fourier fit in black,
# period and epoch in the corner).

const _PALETTE = [:royalblue, :darkorange, :seagreen, :crimson, :purple, :saddlebrown, :deeppink, :gray40, :olive,
                  :darkcyan, :goldenrod, :navy, :firebrick, :teal, :orchid, :black]
_period_text(P, eP) = isnan(eP) ? @sprintf("P = %.4f h (fixed)", 24P) :
                      (d = 24eP < 0.0005 ? 4 : 24eP < 0.005 ? 3 : 2; @sprintf("P = %.*f ± %.*f h", d, 24P, d, 24eP))

"""
    plot(fit::PeriodFit; cycles=1, errorbars=true, residuals=false, title=…)

Phased light curve in the style of the Minor Planet Bulletin: the data of each group (session) in its own colour with
the zero points removed, folded at the fitted period, the Fourier series in black, and the period, its error, the
epoch of phase 0 (JD, light-time corrected if the times were) and the amplitude as an annotation. `cycles = 2` shows
two rotations; `residuals = true` adds a residual panel. Requires Plots.jl.
"""
@recipe function f(fit::PeriodFit; cycles=1, errorbars=true, residuals=false)
    ph = mod.((fit.t .- fit.t0) ./ fit.period, 1.0)
    y = fit.y .- fit.zero_points[fit.g]
    res = y .- fourier_series(fit, ph)
    layout := residuals ? (2, 1) : 1
    size --> (900, residuals ? 620 : 470)
    legend --> :topright
    legendfontsize --> 7
    framestyle --> :box
    for (k, lab) in enumerate(fit.labels)
        m = fit.g .== k
        for c in 0:cycles-1
            @series begin
                seriestype := :scatter
                subplot := 1
                markersize --> 3; markerstrokewidth --> 0.3
                color := _PALETTE[mod1(k, length(_PALETTE))]
                label := c == 0 ? lab : ""
                errorbars && (yerror := fit.σ[m])
                ph[m] .+ c, y[m]
            end
            if residuals
                @series begin
                    seriestype := :scatter; subplot := 2; label := ""; legend := false
                    markersize --> 2.5; markerstrokewidth --> 0
                    color := _PALETTE[mod1(k, length(_PALETTE))]
                    ph[m] .+ c, res[m]
                end
            end
        end
    end
    xs = range(0, cycles; length=600cycles)
    @series begin
        subplot := 1; seriestype := :path; color := :black; linewidth --> 1.8; label := "Fourier order $(fit.order)"
        yflip := true
        lo, hi = extrema(y); rng = max(hi - lo, 0.05)
        ylims --> (lo - 0.42rng, hi + 0.08rng)
        xlims --> (-0.02, cycles + 0.02)
        xlabel --> (residuals ? "" : "Rotational phase")
        ylabel --> "Δmag"
        title --> (isempty(fit.object) ? "Phased light curve" : "Phased light curve: " * fit.object)
        annotations --> [(0.02cycles, lo - 0.40rng, (string(_period_text(fit.period, fit.period_err), "\nJD₀ = ",
                                                   @sprintf("%.5f", fit.t0), "\nAmp = ", @sprintf("%.2f", fit.amplitude),
                                                   " mag"), 9, :left, :top))]
        xs, fourier_series(fit, mod.(xs, 1.0))
    end
    if residuals
        @series begin
            subplot := 2; seriestype := :hline; color := :black; label := ""; yflip := true; legend := false
            xlims --> (-0.02, cycles + 0.02); xlabel --> "Rotational phase"; ylabel --> "O − C"
            [0.0]
        end
    end
end

"""
    plot(blocks::Vector{ALCDEFBlock}; reltime=true)
    plot(block::ALCDEFBlock)

Raw ALCDEF light curves (magnitude against time, one colour per block; time in hours from the start of each session
if `reltime`, else JD). Requires Plots.jl.
"""
@recipe function f(blocks::AbstractVector{ALCDEFBlock}; reltime=true)
    yflip := true; framestyle --> :box; legend --> :outerright
    xlabel --> (reltime ? "hours from session start" : "JD"); ylabel --> "mag"
    for (k, b) in enumerate(blocks)
        @series begin
            seriestype := :scatter; markersize --> 2.5; markerstrokewidth --> 0
            color := _PALETTE[mod1(k, length(_PALETTE))]
            label := string(b["SESSIONDATE", "?"], " ", b["CONTACTNAME", ""])
            yerror := any(isnan, b.err) ? nothing : b.err
            (reltime ? 24 .* (b.jd .- minimum(b.jd)) : b.jd), b.mag
        end
    end
end
@recipe f(b::ALCDEFBlock) = [b]

"""
    periodogramplot(f, power; kind=:gls, periods=Float64[], hours=true)

Periodogram against the period (log axis; `f` in cycles/day, shown in hours if `hours`). `kind` = `:gls` (power),
`:falc` (pass Δχ² = χ²₀ − χ²) or `:pdm` (θ, smaller is better). Vertical lines at the `periods` [h]
(e.g. the best period, its double, the daily alias). Requires Plots.jl.
"""
@userplot PeriodogramPlot
@recipe function f(p::PeriodogramPlot; kind=:gls, periods=Float64[], hours=true)
    f, pw = p.args[1], p.args[2]
    P = hours ? 24 ./ f : 1 ./ f
    xscale --> :log10; framestyle --> :box; legend --> :topright
    xlabel --> (hours ? "period [h]" : "period [d]")
    ylabel --> (kind == :gls ? "GLS power" : kind == :falc ? "FALC Δχ²" : kind == :pdm ? "PDM θ" : "power")
    @series begin
        seriestype := :path
        color --> (kind == :gls ? :royalblue : kind == :falc ? :seagreen : :darkorange); label := ""
        linewidth --> 1
        P, pw
    end
    for (i, m) in enumerate(periods)
        @series begin
            seriestype := :vline; linestyle := :dash; color := :crimson; linewidth := 1
            label := @sprintf("%.3f h", m)
            [hours ? m : m / 24]
        end
    end
end

"""
    phasecurveplot(α, mred; σ=nothing, bands=nothing, ref=nothing, perband=false, fit=nothing)

Phase curve(s): reduced magnitude against the solar phase angle. With `bands`, all bands in one panel shifted by
their offsets relative to `ref` (joint fit by [`fit_phase_curve`](@ref), or `fit`), or one panel per band with its own
H and β (`perband = true`; clipped points as crosses). Requires Plots.jl.
"""
@userplot PhaseCurvePlot
@recipe function f(p::PhaseCurvePlot; σ=nothing, bands=nothing, ref=nothing, perband=false, fit=nothing)
    α, m = p.args[1], p.args[2]
    labels = bands === nothing ? fill("all", length(α)) : collect(bands)
    ub = unique(labels)
    yflip := true; framestyle --> :box
    s = σ === nothing ? ones(length(α)) : σ
    amax = 1.08maximum(α)
    if perband
        layout := length(ub)
        size --> (420min(length(ub), 3), 330cld(length(ub), 3))
        title --> permutedims(["band $b" for b in ub])
        xlabel --> permutedims(fill("phase angle [deg]", length(ub)))
        ylabel --> permutedims(fill("m − 5 log(rΔ)", length(ub)))
        for (j, b) in enumerate(ub)
            k = labels .== b
            pc = fit_phase_curve(α[k], m[k]; σ=σ === nothing ? nothing : s[k])
            @series begin
                subplot := j; seriestype := :scatter; label := "used"; markersize --> 3; markerstrokewidth --> 0.3
                color := :royalblue; yerror := σ === nothing ? nothing : s[k][pc.keep]
                α[k][pc.keep], m[k][pc.keep]
            end
            if any(.!pc.keep)
                @series begin
                    subplot := j; seriestype := :scatter; label := "clipped"; marker := :x; color := :crimson
                    α[k][.!pc.keep], m[k][.!pc.keep]
                end
            end
            @series begin
                subplot := j; seriestype := :path; color := :black; linewidth --> 2; legend --> :bottomright
                label := @sprintf("H = %.2f ± %.2f, β = %.3f ± %.3f", pc.H, pc.H_err, pc.beta, pc.beta_err)
                [0, amax], pc.H .+ pc.beta .* [0, amax]
            end
        end
    else
        pc = fit === nothing ? fit_phase_curve(α, m; σ=σ, bands=bands, ref=ref) : fit
        size --> (900, 520); legend --> :outertopright
        xlabel --> "solar phase angle [deg]"; ylabel --> "reduced magnitude − band offset"
        for (j, b) in enumerate(ub)
            k = (labels .== b) .& pc.keep
            off = get(pc.offsets, b, 0.0)
            @series begin
                seriestype := :scatter; markersize --> 3; markerstrokewidth --> 0.3
                color := _PALETTE[mod1(j, length(_PALETTE))]
                label := bands === nothing ? "data" : @sprintf("%s (%d, %+.2f)", b, count(k), off)
                yerror := σ === nothing ? nothing : s[k]
                α[k], m[k] .- off
            end
        end
        @series begin
            seriestype := :path; color := :black; linewidth --> 2
            label := @sprintf("H_%s = %.2f, β = %.3f", pc.ref, pc.H, pc.beta)
            [0, amax], pc.H .+ pc.beta .* [0, amax]
        end
    end
end

"""
    splithalvesplot(t, y, P; σ=nothing, nbins=10)

Split-halves plot (Harris et al. 2014): the two halves of the light curve folded at `P` overlaid against the phase
within the half-cycle, with binned means and the χ² test of their difference in the title. Requires Plots.jl.
"""
@userplot SplitHalvesPlot
@recipe function f(p::SplitHalvesPlot; σ=nothing, nbins=10)
    t, y, P = p.args[1], p.args[2], p.args[3]
    sh = split_halves(t, y, P; σ=σ, nbins=nbins)
    ph = fold(t, P); h1 = ph .< 0.5; ψ = mod.(ph, 0.5) .* 2
    yflip := true; framestyle --> :box; legend --> :bottomright
    xlabel --> "phase within the half-cycle"; ylabel --> "Δmag"
    title --> @sprintf("split halves at P = %.3f h: χ² = %.1f for %d bins, p = %.2f", 24P, sh.chi2, sh.dof, sh.p)
    for (k, lab, col) in ((h1, "first half", :royalblue), (.!h1, "second half", :darkorange))
        @series begin
            seriestype := :scatter; markersize --> 2.5; markerstrokewidth --> 0; alpha --> 0.5
            color := col; label := lab
            ψ[k], y[k]
        end
    end
    for (mn, er, col) in ((sh.mean1, sh.err1, :royalblue), (sh.mean2, sh.err2, :darkorange))
        @series begin
            seriestype := :path; color := col; linewidth --> 2; marker := :circle; markersize --> 4; label := ""
            yerror := er
            sh.centers, mn
        end
    end
end

"""
    windowplot(t, f; hours=true)

Spectral window of the sampling times `t` against the period (log axis). Requires Plots.jl.
"""
@userplot WindowPlot
@recipe function f(p::WindowPlot; hours=true)
    t, f = p.args[1], p.args[2]
    seriestype := :path; xscale --> :log10; framestyle --> :box; label --> ""; color --> :gray30
    xlabel --> (hours ? "period [h]" : "period [d]"); ylabel --> "spectral window"
    (hours ? 24 ./ f : 1 ./ f), spectral_window(t, f)
end
