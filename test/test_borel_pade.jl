using Test
using Resurgence
using SpecialFunctions: expinti

unsigned_factorials(::Type{T}, N) where {T} = T[factorial(big(k)) for k in 0:N-1]

# Reference value: Stieltjes / Euler's series  S(z) = Σ (-1)^k k! z^k  has
# exact Borel sum  S(z) = (1/z) e^{1/z} E_1(1/z)  for z > 0.
# At z = 1:  S(1) = e · E_1(1) ≈ 0.5963473623231940743…
const STIELTJES_AT_1 = 0.5963473623231940743410784993

stieltjes_coeffs(::Type{T}, N) where {T} = T[(-1)^k * factorial(big(k)) for k in 0:N-1]

@testset "borel_pade — Stieltjes series" begin
    @testset "Float64 reproduces e·E₁(1) at z=1" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 25))   # 25 coeffs avoids Float64 overflow at k≈21
        v = borel_pade(a; n = 10, m = 10, x = 1)
        @test isapprox(v, STIELTJES_AT_1; atol = 1e-8)
    end

    @testset "BigFloat sharpens the answer" begin
        a = stieltjes_coeffs(BigFloat, 41)
        v = borel_pade(a; n = 20, m = 20, x = BigFloat(1))
        # quadgk default tolerances dominate at high precision; this is well
        # below the Float64 precision the same configuration achieves.
        @test abs(v - BigFloat(STIELTJES_AT_1)) < BigFloat("1e-12")
    end

    @testset "return_error returns a tuple" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 25))
        v, err = borel_pade(a; n = 10, m = 10, x = 1, return_error = true)
        @test err ≥ 0
        @test isapprox(v, STIELTJES_AT_1; atol = max(err * 10, 1e-6))
    end

    @testset "regularize_poles only meaningful for x>0" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 25))
        @test_throws ArgumentError borel_pade(a; n = 10, m = 10, x = -1.0,
                                              regularize_poles = true)
    end

    @testset "argument validation" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 5))
        @test_throws ArgumentError borel_pade(a; n = 5, m = 5, x = 1)
    end

    @testset "ComplexF64 coefficient vector" begin
        a = ComplexF64.(stieltjes_coeffs(BigFloat, 25))
        v = borel_pade(a; n = 10, m = 10, x = ComplexF64(1))
        @test v isa Complex
        @test isapprox(v, complex(STIELTJES_AT_1); atol = 1e-8)
    end

    @testset "real coefficients with complex x" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 25))
        v = borel_pade(a; n = 10, m = 10, x = 1.0 + 0.01im)
        @test v isa Complex
        @test isfinite(real(v)) && isfinite(imag(v))
        # small imaginary perturbation ⇒ result close to the real Borel sum
        @test abs(v - STIELTJES_AT_1) < 0.05
    end

    @testset "Complex{BigFloat} propagates" begin
        a = Complex{BigFloat}.(stieltjes_coeffs(BigFloat, 25))
        v = borel_pade(a; n = 10, m = 10, x = Complex{BigFloat}(1))
        @test v isa Complex{BigFloat}
        @test isapprox(v, complex(BigFloat(STIELTJES_AT_1)); atol = BigFloat("1e-8"))
    end

    @testset "regularize_poles requires real positive x" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 25))
        @test_throws ArgumentError borel_pade(a; n = 10, m = 10,
                                              x = 1.0 + 0.01im,
                                              regularize_poles = true)
    end
end

@testset "borel_leroy_pade — Stieltjes series" begin
    @testset "small-b limit recovers Borel sum" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 25))
        # Pure b=0 hits a singular Padé linear system because Borel(Stieltjes)
        # is exactly degree-1 rational; nudge b slightly off zero to lift the
        # degeneracy while keeping the Le Roy weight ≈ 1.
        v = borel_leroy_pade(a; n = 10, m = 10, b = 1e-6, x = 1)
        @test isapprox(v, STIELTJES_AT_1; atol = 1e-5)
    end

    @testset "default b = -1//2 still converges" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 25))
        v = borel_leroy_pade(a; n = 10, m = 10, x = 1)
        # Should still be in the right ballpark — not exact for a generic b.
        @test isapprox(v, STIELTJES_AT_1; atol = 1e-3)
    end

    @testset "ComplexF64 coefficients" begin
        a = ComplexF64.(stieltjes_coeffs(BigFloat, 25))
        v = borel_leroy_pade(a; n = 10, m = 10, b = 1e-6, x = ComplexF64(1))
        @test v isa Complex
        @test isapprox(v, complex(STIELTJES_AT_1); atol = 1e-5)
    end
end

@testset "borel_pade lateral / median / discontinuity" begin
    # Borel-summable case: alternating Stieltjes has its singularity on the
    # negative real axis at t = -1, so no positive-real-axis poles to
    # regularize. Both laterals collapse to the ordinary Borel sum and the
    # discontinuity vanishes.
    @testset "Borel-summable Stieltjes ⇒ disc ≈ 0" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 25))
        ref = borel_pade(a; n = 10, m = 10, x = 1)
        Lp = borel_pade_lateral(a; n = 10, m = 10, x = 1, side = +1)
        Lm = borel_pade_lateral(a; n = 10, m = 10, x = 1, side = -1)
        @test isapprox(Lp, ref; atol = 1e-6)
        @test isapprox(Lm, ref; atol = 1e-6)
        med = borel_pade_median(a; n = 10, m = 10, x = 1)
        disc = borel_pade_discontinuity(a; n = 10, m = 10, x = 1)
        @test isapprox(med, ref; atol = 1e-6)
        @test abs(disc) < 1e-6
    end

    # Non-Borel-summable driver: a_k = k!  ⇒  Borel transform 1/(1-t),
    # pole on positive real axis at t = 1.
    # PV ∫₀^∞ e^{-t}/(1-t) dt = Ei(1)/e   (median).
    # |Discontinuity| = π/e, with sign set by convention below.
    # Here side=+1 shifts the pole to +iε (Plemelj 1/(x-iε) limit), so the
    # lateral L⁺ has imag = -π/e and disc = (L⁺ − L⁻)/(2i) = -π/e.
    @testset "non-Borel-summable a_k = k! ⇒ median = Ei(1)/e, |disc| = π/e" begin
        a = Float64.(unsigned_factorials(BigFloat, 25))
        med = borel_pade_median(a; n = 10, m = 10, x = 1)
        disc = borel_pade_discontinuity(a; n = 10, m = 10, x = 1)
        @test isapprox(real(med), expinti(1) / exp(1); atol = 1e-6)
        @test abs(imag(med)) < 1e-6                       # median is real
        @test isapprox(real(disc), -π / exp(1); atol = 1e-6)
        @test abs(imag(disc)) < 1e-6                       # disc is real
    end

    @testset "lateral(+1) and lateral(-1) are conjugates for real a" begin
        a = Float64.(unsigned_factorials(BigFloat, 25))
        Lp = borel_pade_lateral(a; n = 10, m = 10, x = 1, side = +1)
        Lm = borel_pade_lateral(a; n = 10, m = 10, x = 1, side = -1)
        @test isapprox(Lp, conj(Lm); atol = 1e-6)
    end

    @testset "x > 0 precondition propagates through wrappers" begin
        a = Float64.(unsigned_factorials(BigFloat, 25))
        @test_throws ArgumentError borel_pade_lateral(a; n = 10, m = 10, x = -1.0)
        @test_throws ArgumentError borel_pade_median(a; n = 10, m = 10, x = -1.0)
        @test_throws ArgumentError borel_pade_discontinuity(a; n = 10, m = 10, x = -1.0)
    end

    @testset "Le Roy lateral / median / discontinuity smoke test" begin
        # b ≠ 0 keeps the Padé linear system non-singular; a_k = k! still gives
        # a Borel-Le Roy transform with a positive-real-axis pole, so the
        # lateral wrappers should produce finite complex values whose
        # discontinuity is approximately purely real (= 2 · Im(L⁺)).
        a = Float64.(unsigned_factorials(BigFloat, 25))
        med = borel_leroy_pade_median(a; n = 10, m = 10, b = 1e-3, x = 1)
        disc = borel_leroy_pade_discontinuity(a; n = 10, m = 10, b = 1e-3, x = 1)
        @test isfinite(real(med)) && isfinite(real(disc))
        @test abs(imag(med)) < 1e-3                       # median is ~real
        @test abs(imag(disc)) < 1e-3                       # disc is ~real
        # b → 0 limit should approach the ordinary Borel discontinuity π/e.
        @test isapprox(real(disc), -π / exp(1); atol = 1e-2)
    end
end

@testset "conformal_borel_pade — Stieltjes series" begin
    @testset "Borel singularity at t=-1 ⇒ sing=1 works well" begin
        a = Float64.(stieltjes_coeffs(BigFloat, 25))
        v = conformal_borel_pade(a; n = 10, m = 10, x = 1, sing = 1)
        @test isapprox(v, STIELTJES_AT_1; atol = 1e-6)
    end

    @testset "ComplexF64 coefficients" begin
        a = ComplexF64.(stieltjes_coeffs(BigFloat, 25))
        v = conformal_borel_pade(a; n = 10, m = 10, x = ComplexF64(1), sing = 1)
        @test v isa Complex
        @test isapprox(v, complex(STIELTJES_AT_1); atol = 1e-6)
    end
end
