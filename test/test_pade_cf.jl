using Test
using Resurgence

@testset "pade_cf" begin
    @testset "geometric series 1/(1-z) ↔ [0/1] is exact" begin
        a = ones(Float64, 6)
        v = pade_cf_value(a, 0, 1, 0.5)
        @test v ≈ 1 / (1 - 0.5)
    end

    @testset "exp(z) [3/3] at z = 0.5" begin
        a = [1 / factorial(k) for k in 0:7]
        v = pade_cf_value(a, 3, 3, 0.5)
        @test v ≈ exp(0.5) atol = 1e-6
    end

    @testset "argument validation" begin
        a = ones(5)
        @test_throws ArgumentError pade_cf(a, 3, 3)        # needs length ≥ 7
        @test_throws ArgumentError pade_cf(a, -1, 1)       # n < 0
        @test_throws ArgumentError pade_cf(a, 1, -1)       # m < 0
        @test_throws ArgumentError pade_cf(a, 2, 4)        # asymmetric (qd staircase only)
        @test_throws ArgumentError pade_cf(a, 4, 2)        # asymmetric the other way
        @test_throws ArgumentError pade_cf_value(a, 2, 4, 0.5)
    end

    @testset "BigFloat" begin
        a = [BigFloat(1) / factorial(big(k)) for k in 0:10]
        v = pade_cf_value(a, 4, 4, BigFloat("0.5"))
        @test v isa BigFloat
        # Same accuracy regime as the LU-path Padé[4/4] of exp.
        @test abs(v - exp(BigFloat("0.5"))) < BigFloat("1e-9")
    end

    @testset "Complex eltype" begin
        a = ones(ComplexF64, 6)
        x = 0.5 + 0.1im
        v = pade_cf_value(a, 0, 1, x)
        @test v isa ComplexF64
        @test v ≈ 1 / (1 - x)
    end

    @testset "matches pade() on well-conditioned exp series" begin
        # Both algorithms must produce the same rational function on a series
        # with non-singular Padé matrix. exp at [3/3] is the canonical
        # full-rank diagonal case.
        a = [1 / factorial(k) for k in 0:7]
        for x in (0.1, 0.5, -0.3, 1.2)
            @test pade_cf_value(a, 3, 3, x) ≈ pade_value(a, 3, 3, x)
        end

        # Off-diagonal n + 1 == m case.
        @test pade_cf_value(a, 2, 3, 0.5) ≈ pade_value(a, 2, 3, 0.5)
    end

    @testset "qd division-by-zero is reported, not silently fudged" begin
        # Borel transform of the Stieltjes series is exactly 1/(1+t):
        # b = [(-1)^k] for k = 0..N. The qd ratios q^{(j)}_1 = b[j+2]/b[j+1]
        # are all -1, so e^{(0)}_1 = q^{(1)}_1 - q^{(0)}_1 = 0 — the next
        # rhombus step divides by zero. pade() takes the rank-deficient pinv
        # fallback to recover 2/3; pade_cf has no such fallback and should
        # report the failure cleanly.
        b = Float64[(-1.0)^k for k in 0:6]
        @test_throws ArgumentError pade_cf(b, 2, 2)
        @test_throws ArgumentError pade_cf_value(b, 2, 2, 0.5)
        # And pade() itself still works on this input — the documented split.
        @test pade_value(b, 2, 2, 0.5) ≈ 2 / 3
    end

    @testset "zero leading coefficient is reported" begin
        # a[1] = 0 makes the very first qd ratio undefined.
        a = Float64[0.0, 1.0, 1.0, 1.0]
        @test_throws ArgumentError pade_cf(a, 1, 1)
        @test_throws ArgumentError pade_cf_value(a, 1, 1, 0.5)
    end

    @testset "pade_cf returns Polynomials matching pade()" begin
        a = [1 / factorial(k) for k in 0:7]
        p_cf, q_cf = pade_cf(a, 3, 3)
        p_lu, q_lu = pade(a, 3, 3)
        # The rational function p/q is the same (assuming both are in
        # canonical form with q(0) = 1); compare values rather than
        # individual coefficients to be robust to scalar normalization.
        for x in (0.0, 0.25, 0.7)
            @test p_cf(x) / q_cf(x) ≈ p_lu(x) / q_lu(x)
        end
    end
end
