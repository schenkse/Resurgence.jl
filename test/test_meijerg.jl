using Test
using Resurgence
using SpecialFunctions: expinti

# Reference series: a_k = (-1)^k (k+1)!. Borel transform is exactly 1/(1+t)^2,
# so the Borel sum at z = 1 is
#     ∫_0^∞ e^{-t} / (1+t)^2 dt = 1 - e·E_1(1) ≈ 0.40365263767680594…
# (Stieltjes itself, a_k = (-1)^k k!, has Borel transform 1/(1+t) which is a
# degree-0 rational in our fitting variable — every rational fit is rank-
# deficient, so it isn't a good driver for borel_meijerg. (k+1)! adds one
# nontrivial pole.)
const MEIJERG_REF = 1.0 - exp(1.0) * (-expinti(-1.0))

stieltjes_shifted_coeffs(::Type{T}, N) where {T} =
    T[(-1)^k * factorial(big(k + 1)) for k in 0:N-1]

@testset "borel_meijerg — Stieltjes-shifted series a_k = (-1)^k (k+1)!" begin
    @testset "Float64 odd n hits machine precision" begin
        a = Float64.(stieltjes_shifted_coeffs(BigFloat, 25))
        # Exact rational Borel transform ⇒ smallest l = 1 (n = 3) is enough.
        v = borel_meijerg(a; n = 3, x = 1.0)
        @test v isa Float64
        @test isapprox(v, MEIJERG_REF; atol = 1e-12)
    end

    @testset "Float64 even n triggers normalized branch" begin
        a = Float64.(stieltjes_shifted_coeffs(BigFloat, 25))
        v = borel_meijerg(a; n = 4, x = 1.0)
        @test v isa Float64
        @test isapprox(v, MEIJERG_REF; atol = 1e-10)
    end

    @testset "BigFloat sharpens the answer" begin
        ab = stieltjes_shifted_coeffs(BigFloat, 25)
        vb = borel_meijerg(ab; n = 3, x = BigFloat(1))
        @test vb isa BigFloat
        ref_big = BigFloat(1) - exp(BigFloat(1)) * (-expinti(BigFloat(-1)))
        @test abs(vb - ref_big) < BigFloat("1e-50")
    end

    @testset "ComplexF64 coefficients propagate" begin
        ac = ComplexF64.(stieltjes_shifted_coeffs(BigFloat, 25))
        vc = borel_meijerg(ac; n = 3, x = ComplexF64(1))
        @test vc isa Complex
        @test isapprox(vc, complex(MEIJERG_REF); atol = 1e-12)
    end

    @testset "Real coefficients with complex x give complex result" begin
        a = Float64.(stieltjes_shifted_coeffs(BigFloat, 25))
        v = borel_meijerg(a; n = 3, x = 1.0 + 0.01im)
        @test v isa Complex
        @test isfinite(real(v)) && isfinite(imag(v))
        # Tiny imaginary perturbation ⇒ result close to the real Borel sum.
        @test abs(v - MEIJERG_REF) < 0.01
    end

    @testset "Complex{BigFloat} propagates" begin
        a = Complex{BigFloat}.(stieltjes_shifted_coeffs(BigFloat, 25))
        v = borel_meijerg(a; n = 3, x = Complex{BigFloat}(1))
        @test v isa Complex{BigFloat}
        ref_big = BigFloat(1) - exp(BigFloat(1)) * (-expinti(BigFloat(-1)))
        @test abs(v - complex(ref_big)) < BigFloat("1e-50")
    end

    @testset "argument validation" begin
        a = Float64.(stieltjes_shifted_coeffs(BigFloat, 25))
        @test_throws ArgumentError borel_meijerg(a; n = 2, x = 1.0)
        @test_throws ArgumentError borel_meijerg(a; n = 3, x = 0)
        # length(a) < n + 1
        a_short = Float64.(stieltjes_shifted_coeffs(BigFloat, 3))
        @test_throws ArgumentError borel_meijerg(a_short; n = 3, x = 1.0)
    end

    @testset "agrees with borel_pade on the same series" begin
        a = Float64.(stieltjes_shifted_coeffs(BigFloat, 25))
        v_meijerg = borel_meijerg(a; n = 3, x = 1.0)
        v_pade = borel_pade(a; n = 1, m = 2, x = 1.0)
        # Both methods recover the exact rational ⇒ should agree to ~1e-10.
        @test isapprox(v_meijerg, v_pade; atol = 1e-10)
    end
end
