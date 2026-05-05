using Test
using Resurgence

@testset "shanks" begin
    @testset "Leibniz / Gregory series for π/4" begin
        # 4 ∑ (-1)^k / (2k+1)  → π
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        # plain partial sum is bad; one Shanks application improves it a lot.
        s1 = shanks(a, 5)
        s3 = shanks(a, 5; depth = 3)
        @test abs(s1 - π) < 1e-2
        @test abs(s3 - π) < 1e-5
        @test abs(s3 - π) < abs(s1 - π)
    end

    @testset "argument validation" begin
        a = collect(1.0:5.0)
        @test_throws ArgumentError shanks(a, 1)        # n must be ≥ 2
        @test_throws ArgumentError shanks(a, 5)        # need length ≥ n+1
        @test_throws ArgumentError shanks(a, 2; depth = 0)
    end

    @testset "BigFloat propagates" begin
        a = BigFloat[(-1)^k * BigFloat(1) / (2k + 1) for k in 0:20]
        s = shanks(4 * a, 5; depth = 2)
        @test s isa BigFloat
        @test abs(s - BigFloat(π)) < 1e-3
    end

    @testset "converged tail does not divide by zero" begin
        # Series with all-zero tail: partial-sum differences are identically 0,
        # the Shanks denominator vanishes. Expect the partial sum (1.0), not NaN.
        a = [1.0, 0.0, 0.0, 0.0]
        s1 = shanks(a, 2)
        @test isfinite(s1)
        @test s1 ≈ 1.0

        # Same with iterated depth (needs n ≥ depth+1 for the iter buffer).
        a2 = zeros(8); a2[1] = 1.0
        s2 = shanks(a2, 3; depth = 2)
        @test isfinite(s2)
        @test s2 ≈ 1.0
    end

    @testset "linear partial sums (constant tail) return the partial sum" begin
        # Constant a[k] = c ⇒ Aₙ linear in n ⇒ d1 == d2; Shanks denominator
        # vanishes. Implementation falls back to the partial sum.
        a = fill(0.5, 6)
        s = shanks(a, 3)
        @test isfinite(s)
        @test s ≈ sum(a[1:4])
    end

    @testset "Complex eltype" begin
        # Same Leibniz series, scaled by a complex constant; Shanks is
        # entirely linear in `a`, so the result must be (1 + 0.1im) × real.
        N = 30
        a_real = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        a_cplx = ComplexF64.(a_real) .* (1 + 0.1im)
        s_real = shanks(a_real, 5; depth = 3)
        s_cplx = shanks(a_cplx, 5; depth = 3)
        @test s_cplx isa ComplexF64
        @test s_cplx ≈ s_real * (1 + 0.1im)
    end
end
