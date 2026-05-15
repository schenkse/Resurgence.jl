using Test
using Resurgence

@testset "theta_brezinski" begin
    @testset "Leibniz / Gregory series for π" begin
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        v1 = theta_brezinski(a, 5; depth = 1)
        v3 = theta_brezinski(a, 5; depth = 3)
        @test abs(v1 - π) < 1e-2
        @test abs(v3 - π) < 1e-5
        @test abs(v3 - π) < abs(v1 - π)
    end

    @testset "argument validation" begin
        a = collect(1.0:6.0)
        @test_throws ArgumentError theta_brezinski(a, 1; depth = 0)
        @test_throws ArgumentError theta_brezinski(a, 0; depth = 1)
        @test_throws ArgumentError theta_brezinski(a, 1; depth = 2)  # need n+3·depth ≤ length
    end

    @testset "BigFloat propagates" begin
        N = 30
        a = BigFloat[4 * (-1)^k / (2k + 1) for k in 0:N-1]
        v = theta_brezinski(a, 5; depth = 3)
        @test v isa BigFloat
        @test abs(v - BigFloat(π)) < 1e-4
    end

    @testset "Complex linearity" begin
        N = 30
        a_real = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        a_cplx = ComplexF64.(a_real) .* (1 + 0.1im)
        v_real = theta_brezinski(a_real, 5; depth = 3)
        v_cplx = theta_brezinski(a_cplx, 5; depth = 3)
        @test v_cplx isa ComplexF64
        @test v_cplx ≈ v_real * (1 + 0.1im)
    end

    @testset "converged tail does not divide by zero" begin
        a = zeros(Float64, 12); a[1] = 1.0
        v = theta_brezinski(a, 2; depth = 2)
        @test isfinite(v)
        @test v ≈ 1.0
    end
end
