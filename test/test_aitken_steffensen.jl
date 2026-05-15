using Test
using Resurgence

@testset "aitken_steffensen" begin
    @testset "Dottie number — fixed point of cos" begin
        # cos has a unique real fixed point at the Dottie number ≈ 0.7390851332151607.
        dottie = 0.7390851332151607
        x = aitken_steffensen(cos, 0.5)
        @test isfinite(x)
        @test abs(x - dottie) < 1e-12
        # Reached cosine fixed point: x ≈ cos(x).
        @test abs(x - cos(x)) < 1e-12
    end

    @testset "quadratic g(x) = (x² + 2)/3 has fixed points 1 and 2" begin
        g(x) = (x^2 + 2) / 3
        # Starting near 1 converges to 1 (the attractive fixed point).
        x = aitken_steffensen(g, 0.9)
        @test abs(x - 1.0) < 1e-12
        @test abs(g(x) - x) < 1e-12
    end

    @testset "BigFloat precision" begin
        # Dottie number to high precision: solve x = cos(x) in BigFloat.
        x = aitken_steffensen(cos, BigFloat("0.5"))
        @test x isa BigFloat
        @test abs(x - cos(x)) < BigFloat(10)^(-50)
    end

    @testset "non-convergence raises" begin
        # g(x) = x + 1 has no fixed point; iteration diverges.
        @test_throws ArgumentError aitken_steffensen(x -> x + 1, 0.0;
                                                     maxiter = 5)
        # Validation: maxiter must be ≥ 1.
        @test_throws ArgumentError aitken_steffensen(cos, 0.5; maxiter = 0)
    end

    @testset "already at the fixed point returns immediately" begin
        # cos's fixed point: starting exactly there should be detected by the
        # denominator check.
        dottie = 0.7390851332151607
        x = aitken_steffensen(cos, dottie)
        @test abs(x - dottie) < 1e-12
    end
end
