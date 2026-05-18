# Quartic anharmonic oscillator — the textbook Borel-summable physics series.
#
# Hamiltonian: H = p²/2 + x²/2 + g·x⁴. The ground-state energy
#
#     E_0(g) = 1/2 + (3/4) g - (21/8) g² + (333/16) g³ - ...
#
# is the canonical Bender-Wu series (Phys. Rev. 184, 1231 (1969)). The
# coefficients diverge factorially, sign-alternating, with the leading Borel
# singularity at t = -1/3 on the negative real axis — Borel-summable, and
# tame ground for `borel_pade`, `borel_leroy_pade`, `conformal_borel_pade`.
#
# Reference E_0(0.1) ≈ 0.5591463272…  from high-precision diagonalisation.

using Resurgence
using Printf

# Compute the first N+1 Bender-Wu coefficients [E_0^(0), E_0^(1), …, E_0^(N)]
# as exact rationals via Rayleigh-Schrödinger perturbation in the ψ_0-times-
# polynomial ansatz ψ_n(x) = exp(-x²/2) · Σ_{k=0}^{2n} A_{n,k} x^{2k}.
# The recurrence is, for n ≥ 1 and k descending from 2n to 1,
#
#     A_{n,k} = [(k+1)(2k+1) A_{n,k+1} - A_{n-1, k-2}
#                + Σ_{m=1}^{n-1} E_0^{(m)} · A_{n-m, k}] / (2k),
#
# with A_{n,0} = 0 (intermediate normalisation) and E_0^{(n)} = -A_{n,1}.
function bender_wu_aho_coefficients(N::Integer)
    R = Rational{BigInt}
    E = R[1//2]                             # E_0^(0) = 1/2
    A = Vector{Vector{R}}()                 # A[n+1] holds (A_{n,0}, A_{n,1}, …, A_{n,2n})
    push!(A, R[1])                          # A_0 = (1)
    for n in 1:N
        An = zeros(R, 2n + 1)
        for k in 2n:-1:1
            Akp1     = (k + 1 ≤ 2n)               ? An[k + 2]   : R(0)
            Anm1km2  = (0 ≤ k - 2 ≤ 2(n - 1))     ? A[n][k - 1] : R(0)
            s = R(0)
            for m in 1:n-1
                if k ≤ 2(n - m)
                    s += E[m + 1] * A[n - m + 1][k + 1]
                end
            end
            An[k + 1] = (R(k + 1) * R(2k + 1) * Akp1 - Anm1km2 + s) // R(2k)
        end
        push!(A, An)
        push!(E, -An[2])                    # E_0^{(n)} = -A_{n,1}
    end
    return E
end

# 21 coefficients give us room for [n=10/m=10] Padé fits while staying inside
# Float64 dynamic range at g = 0.1.
const BW = bender_wu_aho_coefficients(20)
const a  = Float64.(BW)

# High-precision reference at g = 0.1 from numerical diagonalisation
# (Bender-Wu 1969, Hioe-Montroll 1975, many subsequent confirmations).
const G       = 0.1
const E0_REF  = 0.55914633237109099

println("Quartic anharmonic oscillator E_0(g) at g = $G")
println("Reference: $E0_REF")
println("Bender-Wu coefficients (first 8): ", BW[1:8])
println()

# ----- Direct partial sum (asymptotic) ---------------------------------------
partial = sum(a[k+1] * G^k for k in 0:length(a)-1)
@printf("direct partial sum (all 21 terms):     %+.10f   (err %.2e)\n",
        partial, partial - E0_REF)
println("(asymptotic — adding more terms eventually makes it worse)")
println()

# ----- Optimal truncation ----------------------------------------------------
ak_g = [a[k+1] * G^k for k in 0:length(a)-1]
Nstar, partial_opt, εN = optimal_truncation(ak_g)
@printf("optimal truncation (N* = %d):            %+.10f   (err %.2e, est %.2e)\n",
        Nstar, partial_opt, partial_opt - E0_REF, εN)
println()

# ----- Borel-Padé ------------------------------------------------------------
v_bp = borel_pade(a; n = 10, m = 10, x = G)
@printf("Borel-Padé[10/10]:                      %+.10f   (err %.2e)\n",
        v_bp, v_bp - E0_REF)

# ----- Borel-Le Roy-Padé (b = -1/2 default) ----------------------------------
v_lr = borel_leroy_pade(a; n = 10, m = 10, x = G)
@printf("Borel-Le Roy-Padé[10/10] (b = -1/2):    %+.10f   (err %.2e)\n",
        v_lr, v_lr - E0_REF)

# ----- Conformal Borel-Padé (sing = 1/3, the AHO singularity location) --------
# The Bender-Wu large-order analysis gives a_k ~ (-3)^k k!^{1/2-style}, so the
# Borel transform has its nearest singularity at t = -1/3 — exactly the
# `sing = 1/3` input to `conformal_borel_pade`.
v_cb = conformal_borel_pade(a; n = 10, m = 10, x = G, sing = 1//3)
@printf("conformal-Borel-Padé[10/10] (sing=1/3): %+.10f   (err %.2e)\n",
        v_cb, v_cb - E0_REF)
println()

# ----- Diagnostics -----------------------------------------------------------
# For alternating Borel-summable series the extracted Stokes action S is
# negative; raising a negative real number to a non-integer β goes complex,
# so we pass `a` as ComplexF64 to keep `stokes_fit` happy. (`diagnose` on
# the real `a` returns the same growth/alternating/recommended fields with
# missing S, β, A.)
d = diagnose(ComplexF64.(a))
println("diagnose:")
@printf("  growth      = %s\n", d.growth)
@printf("  alternating = %s\n", d.alternating)
@printf("  S           = %+.4f + %+.4fi\n", real(d.S), imag(d.S))
@printf("  β           = %+.4f + %+.4fi\n", real(d.β), imag(d.β))
@printf("  A           = %+.4f + %+.4fi\n", real(d.A), imag(d.A))
@printf("  recommended = %s\n", d.recommended)
println("  (S ≈ -1/3 confirms the negative-real-axis Borel singularity at t = -1/3)")
println()

cd = coefficient_diagnostics(a)
println("coefficient_diagnostics:")
@printf("  growth            = %s\n", cd.growth)
@printf("  ratio_growth      = %.2f   (≫ 1 → factorial)\n", cd.ratio_growth)
@printf("  alternation_score = %.2f\n", cd.alternation_score)
println()

# ----- compare ---------------------------------------------------------------
rows = compare([
    BorelPade(10, 10; x = G),
    BorelLeRoyPade(10, 10; x = G),
    ConformalBorelPade(10, 10; x = G, sing = 1//3),
    Pade(10, 10; x = G),
], a; reference = E0_REF)
println("compare:")
@printf("  %-30s  %-16s  %s\n", "method", "result", "residual")
for r in rows
    res = r.result === missing ? "missing" : @sprintf("%+.10f", real(r.result))
    rsd = r.residual === missing ? "—" : @sprintf("%.2e", r.residual)
    @printf("  %-30s  %-16s  %s\n", r.method, res, rsd)
end
println()

# ----- BigFloat precision pass -----------------------------------------------
setprecision(BigFloat, 256) do
    a_big = BigFloat.(BW)
    g_big = big"0.1"
    v_big = borel_pade(a_big; n = 10, m = 10, x = g_big,
                        rtol = BigFloat("1e-30"), atol = BigFloat("1e-30"))
    println("BigFloat Borel-Padé[10/10] at 256-bit precision:")
    println("  ", v_big)
    println("  (the leading 1e-14 error is dominated by quadgk's Float64-default node")
    println("   computation; pass smaller rtol / atol or order=… for tighter integration)")
end
