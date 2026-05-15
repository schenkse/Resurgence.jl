using Test
using Resurgence

@testset "wynn_eps" begin
    @testset "Leibniz / Gregory series for π/4" begin
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        s1 = wynn_eps(a, 5)
        s3 = wynn_eps(a, 5; depth = 3)
        @test abs(s1 - π) < 1e-2
        @test abs(s3 - π) < 1e-5
        @test abs(s3 - π) < abs(s1 - π)
    end

    @testset "depth=1 matches single Aitken-Δ²" begin
        # Column 2 of the ε-tableau coincides with one Shanks step.
        # Higher depths differ entry-by-entry (ε realises k-th-order Shanks
        # via the determinantal form; `shanks` iterates the Δ² operator),
        # so only depth=1 is an algebraic identity.
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        @test wynn_eps(a, 5; depth = 1) ≈ shanks(a, 5; depth = 1)
    end

    @testset "exact on a single geometric tail" begin
        # S_n = L + c r^n with r = 1/2, L = 2, c = -1.
        # Column 2 of the ε-tableau is exact (Δ² annihilates one geometric tail).
        r, L, c = 0.5, 2.0, -1.0
        N = 6
        S = Float64[L + c * r^n for n in 1:N]
        a = vcat(S[1], diff(S))
        @test wynn_eps(a, 4; depth = 1) ≈ L
    end

    @testset "argument validation" begin
        a = collect(1.0:5.0)
        @test_throws ArgumentError wynn_eps(a, 1)              # n must be ≥ 2
        @test_throws ArgumentError wynn_eps(a, 5)              # need length ≥ n+depth
        @test_throws ArgumentError wynn_eps(a, 2; depth = 0)
        @test_throws ArgumentError wynn_eps(a, 2; depth = 2)   # n < depth+1 → underfilled column
    end

    @testset "BigFloat propagates" begin
        a = BigFloat[(-1)^k * BigFloat(1) / (2k + 1) for k in 0:20]
        s = wynn_eps(4 * a, 5; depth = 2)
        @test s isa BigFloat
        @test abs(s - BigFloat(π)) < 1e-3
    end

    @testset "Complex eltype" begin
        N = 30
        a_real = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        a_cplx = ComplexF64.(a_real) .* (1 + 0.1im)
        s_real = wynn_eps(a_real, 5; depth = 3)
        s_cplx = wynn_eps(a_cplx, 5; depth = 3)
        @test s_cplx isa ComplexF64
        @test s_cplx ≈ s_real * (1 + 0.1im)
    end

    @testset "converged tail does not divide by zero" begin
        # All-zero tail: consecutive partial sums coincide, the ε denominator
        # vanishes. Expect the partial sum, not NaN.
        a = [1.0, 0.0, 0.0, 0.0]
        s1 = wynn_eps(a, 2)
        @test isfinite(s1)
        @test s1 ≈ 1.0

        a2 = zeros(8); a2[1] = 1.0
        s2 = wynn_eps(a2, 3; depth = 2)
        @test isfinite(s2)
        @test s2 ≈ 1.0
    end
end
