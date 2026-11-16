# Small statistical helpers (no SpecialFunctions dependency).

# log Γ(x) for x > 0 (Lanczos approximation, g = 7, n = 9; relative error ~1e-15)
const _LANCZOS = (0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
                  -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
                  1.5056327351493116e-7)
function _loggamma(x::Real)
    x < 0.5 && return log(π / abs(sin(π * x))) - _loggamma(1 - x)
    x -= 1
    a = _LANCZOS[1]
    t = x + 7.5
    for i in 1:8
        a += _LANCZOS[i+1] / (x + i)
    end
    return 0.5 * log(2π) + (x + 0.5) * log(t) - t + log(a)
end

# regularised lower incomplete gamma P(a, x): series for x < a + 1, continued fraction otherwise
function _gammap(a::Real, x::Real)
    x <= 0 && return 0.0
    if x < a + 1
        ap = a; s = 1 / a; d = s
        for _ in 1:1000
            ap += 1; d *= x / ap; s += d
            abs(d) < abs(s) * 1e-15 && break
        end
        return s * exp(-x + a * log(x) - _loggamma(a))
    else
        b = x + 1 - a; c = 1 / 1e-300; d = 1 / b; h = d
        for i in 1:1000
            an = -i * (i - a); b += 2
            d = an * d + b; abs(d) < 1e-300 && (d = 1e-300)
            c = b + an / c; abs(c) < 1e-300 && (c = 1e-300)
            d = 1 / d; del = d * c; h *= del
            abs(del - 1) < 1e-15 && break
        end
        return 1 - exp(-x + a * log(x) - _loggamma(a)) * h
    end
end

"""
    chi2_cdf(x, k) -> Float64

Cumulative distribution function of the χ² distribution with `k` degrees of freedom; the p-value of a χ² statistic is
`1 - chi2_cdf(x, k)`.
"""
chi2_cdf(x::Real, k::Real) = (k > 0 || throw(ArgumentError("k must be positive")); _gammap(k / 2, x / 2))

# weighted linear least squares: coefficients, covariance (unscaled), χ²
function _wls(X::AbstractMatrix, y::AbstractVector, w::AbstractVector)
    sw = sqrt.(w)
    A = sw .* X
    F = qr(A)
    c = F \ (sw .* y)
    R = F.R
    C = inv(R' * R)
    χ2 = sum(abs2, sw .* y .- A * c)
    return c, C, χ2
end

_weights(σ, n) = σ === nothing ? ones(n) : (length(σ) == n || throw(DimensionMismatch("σ has the wrong length"));
                                            all(>(0), σ) || throw(ArgumentError("σ must be positive")); 1 ./ σ .^ 2)
