using Test
using Resurgence
using SpecialFunctions: expinti

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
