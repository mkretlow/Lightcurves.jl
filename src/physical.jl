# Physical quantities from photometry.

const _K_DIAM = 1329.0      # km: D = 1329 km / √p_V × 10^(−H_V/5) (e.g. Harris & Lagerros 2002)

"""
    geometric_albedo(H, D; H_err=0.0, D_err=0.0) -> (p, p_err)

Geometric albedo from the absolute magnitude `H` (in the band of the albedo, usually V) and the diameter `D` [km]:
p = (1329 km / D)² 10^(−0.4 H), with first-order error propagation.
"""
function geometric_albedo(H::Real, D::Real; H_err::Real=0.0, D_err::Real=0.0)
    D > 0 || throw(ArgumentError("D must be positive"))
    p = (_K_DIAM / D)^2 * 10^(-0.4H)
    return p, p * sqrt((2D_err / D)^2 + (0.4log(10) * H_err)^2)
end

"""
    diameter(H, p; H_err=0.0, p_err=0.0) -> (D, D_err)

Diameter [km] from `H` and the geometric albedo `p`: D = 1329 km / √p × 10^(−H/5).
"""
function diameter(H::Real, p::Real; H_err::Real=0.0, p_err::Real=0.0)
    p > 0 || throw(ArgumentError("p must be positive"))
    D = _K_DIAM / sqrt(p) * 10^(-H / 5)
    return D, D * sqrt((p_err / (2p))^2 + (log(10) / 5 * H_err)^2)
end

"""
    ellipsoid_amplitude(a, b, c, θ) -> Δm [mag]

Peak-to-peak light-curve amplitude of a triaxial ellipsoid with semi-axes a ≥ b ≥ c (rotation about c), seen at the
aspect angle θ [rad] between the spin axis and the line of sight, at zero phase angle (geometric scattering):
Δm = 2.5 log10(a/b) − 1.25 log10[(a² cos²θ + c² sin²θ) / (b² cos²θ + c² sin²θ)]
(e.g. Binzel et al. 1989, Asteroids II). Equator-on (θ = π/2): 2.5 log10(a/b); pole-on (θ = 0): 0.
"""
function ellipsoid_amplitude(a::Real, b::Real, c::Real, θ::Real)
    (a >= b >= c > 0) || throw(ArgumentError("need a ≥ b ≥ c > 0"))
    c2, s2 = cos(θ)^2, sin(θ)^2
    return 2.5log10(a / b) - 1.25log10((a^2 * c2 + c^2 * s2) / (b^2 * c2 + c^2 * s2))
end
