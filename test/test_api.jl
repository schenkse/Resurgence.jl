using Test
using Resurgence

@testset "unified resum API" begin
    a_stieltjes = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
    direct = borel_pade(a_stieltjes; n = 10, m = 10, x = 1)

    @testset "BorelPade dispatches to borel_pade" begin
        v = resum(BorelPade(10, 10), a_stieltjes)
        @test v ≈ direct
    end

    @testset "BorelPadeLateral / BorelPadeMedian dispatch" begin
        a_unsigned = Float64[Float64(factorial(big(k))) for k in 0:24]
        v_lat_tag = resum(BorelPadeLateral(10, 10; x = 1, side = +1), a_unsigned)
        v_lat_direct = borel_pade_lateral(a_unsigned; n = 10, m = 10, x = 1, side = +1)
        @test v_lat_tag ≈ v_lat_direct
        v_med_tag = resum(BorelPadeMedian(10, 10; x = 1), a_unsigned)
        v_med_direct = borel_pade_median(a_unsigned; n = 10, m = 10, x = 1)
        @test v_med_tag ≈ v_med_direct
    end

    @testset "BorelLeRoyPade dispatches to borel_leroy_pade" begin
        # b ≠ 0 keeps the Padé linear system non-singular on this exact-rational
        # Borel transform.
        v = resum(BorelLeRoyPade(10, 10; b = -0.25), a_stieltjes)
        v_direct = borel_leroy_pade(a_stieltjes; n = 10, m = 10, b = -0.25, x = 1)
        @test v ≈ v_direct
    end

    @testset "ConformalBorelPade dispatches" begin
        v = resum(ConformalBorelPade(10, 10; sing = 1.0), a_stieltjes)
        v_direct = conformal_borel_pade(a_stieltjes; n = 10, m = 10, x = 1, sing = 1.0)
        @test v ≈ v_direct
    end

    @testset "HermitePade dispatches to hermite_pade_value (with and without branch)" begin
        # exp(z/2) Taylor coefficients: non-degenerate driver for hermite_pade
        # (f² = exp(z) has no polynomial relation in f).
        N = 11
        c = Vector{Float64}(undef, N)
        c[1] = 1.0
        for k in 1:N-1
            c[k+1] = c[k] / (2k)
        end
        z = 0.3
        v_default = resum(HermitePade(3, 3, 3; x = z), c)
        @test v_default ≈ hermite_pade_value(c, 3, 3, 3, z)
        v_plus  = resum(HermitePade(3, 3, 3; x = z, branch = +1), c)
        v_minus = resum(HermitePade(3, 3, 3; x = z, branch = -1), c)
        @test v_plus  ≈ hermite_pade_value(c, 3, 3, 3, z; branch = +1)
        @test v_minus ≈ hermite_pade_value(c, 3, 3, 3, z; branch = -1)
    end

    @testset "ConformalBorelPadePair dispatches to conformal_borel_pade_pair" begin
        # B(t) = 1/(1+t²) driver with ±i Borel singularities.
        a = Float64[isodd(k) ? 0.0 : (-1.0)^(k ÷ 2) * Float64(factorial(big(k))) for k in 0:24]
        v_tag = resum(ConformalBorelPadePair(10, 10; sing = 1.0), a)
        v_direct = conformal_borel_pade_pair(a; n = 10, m = 10, x = 1, sing = 1.0)
        @test v_tag ≈ v_direct
    end

    @testset "MeijerG dispatches to borel_meijerg" begin
        # Stieltjes itself is rank-degenerate for borel_meijerg; use the
        # shifted (k+1)! variant to drive a clean rational fit.
        a_shifted = Float64[(-1.0)^k * Float64(factorial(big(k + 1))) for k in 0:24]
        v_tag = resum(MeijerG(3), a_shifted)
        v_direct = borel_meijerg(a_shifted; n = 3, x = 1)
        @test v_tag ≈ v_direct
    end

    @testset "Pade dispatches to pade_value" begin
        a = ones(Float64, 6)
        v = resum(Pade(0, 1; x = 0.5), a)
        @test v ≈ 1 / (1 - 0.5)
    end

    @testset "Pade and BorelPade share (n, m) convention" begin
        # Asymmetric (n, m) — assert positional ordering is observable and that
        # both APIs interpret the first arg as the numerator degree.
        @test pade_value(a_stieltjes, 8, 12, 0.5) != pade_value(a_stieltjes, 12, 8, 0.5)
        @test resum(Pade(8, 12; x = 0.5), a_stieltjes) ≈ pade_value(a_stieltjes, 8, 12, 0.5)
        @test resum(BorelPade(8, 12), a_stieltjes) ≈ borel_pade(a_stieltjes; n = 8, m = 12)
    end

    @testset "Shanks dispatches with depth" begin
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        @test resum(Shanks(5; depth = 3), a) ≈ shanks(a, 5; depth = 3)
    end

    @testset "WynnEps dispatches with depth" begin
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        @test resum(WynnEps(5; depth = 3), a) ≈ wynn_eps(a, 5; depth = 3)
    end

    @testset "BrezinskiTheta dispatches with depth" begin
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        @test resum(BrezinskiTheta(5; depth = 3), a) ≈
              theta_brezinski(a, 5; depth = 3)
    end

    @testset "BrezinskiRho dispatches with depth" begin
        a = Float64[1 / k^2 for k in 1:30]
        @test resum(BrezinskiRho(10; depth = 3), a) ≈
              rho_brezinski(a, 10; depth = 3)
    end

    @testset "Cesaro dispatches with depth" begin
        N = 100
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        @test resum(Cesaro(N; depth = 2), a) ≈ cesaro(a, N; depth = 2)
    end

    @testset "Abel dispatches with x" begin
        a = ones(Float64, 8)
        @test resum(Abel(; x = 0.5), a) ≈ abel(a; x = 0.5)
    end

    @testset "Richardson dispatches with depth" begin
        a = Float64[1 / k^2 for k in 1:30]
        @test resum(Richardson(10; depth = 3), a) ≈ richardson(a, 10; depth = 3)
    end

    @testset "Levin dispatches with depth/variant/β" begin
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        @test resum(Levin(3; depth = 20, variant = :u), a) ≈
              levin(a, 3; depth = 20, variant = :u)
        @test resum(Levin(3; depth = 20, variant = :t), a) ≈
              levin(a, 3; depth = 20, variant = :t)
        @test resum(Levin(3; variant = :u), a) ≈ levin(a, 3; variant = :u)
    end

    @testset "Weniger dispatches with depth/β" begin
        a = Float64[(-1.0)^k * Float64(factorial(big(k))) for k in 0:24]
        @test resum(Weniger(3; depth = 15), a) ≈ weniger(a, 3; depth = 15)
        @test resum(Weniger(3), a) ≈ weniger(a, 3)
    end

    @testset "SidiS dispatches with depth/variant/β" begin
        N = 30
        a = Float64[4 * (-1.0)^k / (2k + 1) for k in 0:N-1]
        @test resum(SidiS(3; depth = 15, variant = :u), a) ≈
              sidi_s(a, 3; depth = 15, variant = :u)
        @test resum(SidiS(3), a) ≈ sidi_s(a, 3)
    end

    @testset "kwargs forwarded to underlying call" begin
        # rtol forwarded; loose tol shouldn't make this fail catastrophically.
        v = resum(BorelPade(10, 10; rtol = 1e-3), a_stieltjes)
        @test isapprox(v, 0.5963473623; atol = 1e-3)
    end
end
