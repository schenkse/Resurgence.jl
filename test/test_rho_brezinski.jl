using Test
using Resurgence

@testset "rho_brezinski" begin
    @testset "harmonic-square Σ 1/k² → π²/6" begin
        N = 30
        a = Float64[1 / k^2 for k in 1:N]
        v1 = rho_brezinski(a, 10; depth = 1)
        v3 = rho_brezinski(a, 10; depth = 3)
        @test abs(v1 - π^2 / 6) < 1e-4
        @test abs(v3 - π^2 / 6) < 1e-9
        @test abs(v3 - π^2 / 6) < abs(v1 - π^2 / 6)
    end

    @testset "logarithmically convergent Σ 1/k^1.5 — accelerates plain sum" begin
        # Slow non-integer power-law tail; ρ accelerates but does not reach
        # machine precision — the meaningful claim is "much better than the
        # corresponding partial sum on the same data window".
        N = 40
        ref = 2.6123753486854883  # ζ(1.5)
        a = Float64[1 / k^1.5 for k in 1:N]
        d = 5
        v = rho_brezinski(a, 5; depth = d)
        plain = sum(a[1:5+2*d])
        @test abs(v - ref) < abs(plain - ref) / 5
    end

    @testset "argument validation" begin
        a = collect(1.0:6.0)
        @test_throws ArgumentError rho_brezinski(a, 1; depth = 0)
        @test_throws ArgumentError rho_brezinski(a, 0; depth = 1)
        @test_throws ArgumentError rho_brezinski(a, 3; depth = 2)  # need n+2·depth ≤ length
    end

    @testset "BigFloat propagates" begin
        N = 30
        a = BigFloat[1 / BigFloat(k)^2 for k in 1:N]
        v = rho_brezinski(a, 10; depth = 3)
        @test v isa BigFloat
        @test abs(v - BigFloat(π)^2 / 6) < 1e-9
    end

    @testset "Complex linearity" begin
        N = 30
        a_real = Float64[1 / k^2 for k in 1:N]
        a_cplx = ComplexF64.(a_real) .* (1 + 0.1im)
        v_real = rho_brezinski(a_real, 10; depth = 3)
        v_cplx = rho_brezinski(a_cplx, 10; depth = 3)
        @test v_cplx isa ComplexF64
        @test v_cplx ≈ v_real * (1 + 0.1im)
    end

    @testset "converged tail does not divide by zero" begin
        a = zeros(Float64, 10); a[1] = 1.0
        v = rho_brezinski(a, 2; depth = 2)
        @test isfinite(v)
        @test v ≈ 1.0
    end
end
