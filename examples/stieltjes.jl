# Stieltjes / Euler's series:  S(z) = Σ_{k≥0} (-1)^k k! z^k
#
# Divergent for every z ≠ 0, but Borel summable to
#     S(z) = (1/z) e^{1/z} E₁(1/z),    z > 0.
# At z = 1 this evaluates to  e · E₁(1) ≈ 0.5963473623…
#
# Run with
#     julia --project=. examples/stieltjes.jl
# from the repository root.

using Resurgence
using Printf

const STIELTJES_AT_1 = 0.5963473623231940743410784993

println("Stieltjes / Euler's series at z = 1")
println("Reference value: ", STIELTJES_AT_1)
println()

# 25 coefficients in Float64 (factorials beyond k≈20 lose precision but stay
# representable up to k≈170 in Float64).
a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]

# Direct partial sum (asymptotic — gets worse with more terms).
partial = sum(a)

# Optimal-truncation
Nstar, partial_opt, εN = optimal_truncation(a)

# Various resummation methods.
v_borel_pade = borel_pade(a; n = 10, m = 10, x = 1)
v_borel_lr   = borel_leroy_pade(a; n = 10, m = 10, x = 1)             # b = -1//2 by default
v_conformal  = conformal_borel_pade(a; n = 10, m = 10, x = 1, sing = 1)

@printf "  partial sum (all 25 terms)         : %+.6e   (err %.2e)\n" partial (partial - STIELTJES_AT_1)
@printf "  optimal truncation (N* = %2d)        : %+.6e   (err %.2e, est %.2e)\n" Nstar partial_opt (partial_opt - STIELTJES_AT_1) εN
@printf "  Borel–Padé[10/10]                  : %+.6e   (err %.2e)\n" v_borel_pade (v_borel_pade - STIELTJES_AT_1)
@printf "  Borel–Le Roy–Padé[10/10] (b=-1/2)  : %+.6e   (err %.2e)\n" v_borel_lr (v_borel_lr - STIELTJES_AT_1)
@printf "  conformal-Borel–Padé[10/10]        : %+.6e   (err %.2e)\n" v_conformal (v_conformal - STIELTJES_AT_1)

println()
println("Same problem at BigFloat precision (41 coefficients, [20/20]):")
ab = BigFloat[(-1)^k * factorial(big(k)) for k in 0:40]
vb = borel_pade(ab; n = 20, m = 20, x = BigFloat(1))
@printf "  Borel–Padé[20/20] (BigFloat)       : %s\n" string(vb)
@printf "  error vs reference                 : %s\n" string(abs(vb - BigFloat(STIELTJES_AT_1)))
