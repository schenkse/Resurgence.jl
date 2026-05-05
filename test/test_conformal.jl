using Test
using Resurgence

@testset "conformal" begin
    @testset "map / inverse round-trip" begin
        for t in (0.1, 0.5, 1.0, 5.0, 100.0)
            w = conformal_map(t; a = 1)
            @test inverse_conformal_map(w; a = 1) ≈ t
        end
    end

    @testset "map at t = 0 returns 0" begin
        @test conformal_map(0.0; a = 1) ≈ 0.0
    end

    @testset "argument validation" begin
        @test_throws ArgumentError conformal_map(1.0; a = -1.0)
        @test_throws ArgumentError inverse_conformal_map(0.5; a = 0.0)
    end

    @testset "real t in the cut region returns complex" begin
        # t < -a sits on the cut; the map should return a complex value, not throw.
        w = conformal_map(-2.0; a = 1.0)
        @test w isa Complex
        @test inverse_conformal_map(w; a = 1.0) ≈ -2.0
    end

    @testset "real t in the regular region stays real" begin
        # t > -a: the result should still be real-valued (no spurious promotion).
        w = conformal_map(0.5; a = 1.0)
        @test w isa Real
    end

    @testset "complex t passes through" begin
        w = conformal_map(0.5 + 0.1im; a = 1.0)
        @test w isa Complex
        @test inverse_conformal_map(w; a = 1.0) ≈ 0.5 + 0.1im
    end

    @testset "conformal_reseries reproduces series at small w" begin
        # B(t) = 1/(1+t)  has coefficients (-1)^k.
        N = 12
        b = Float64[(-1.0)^k for k in 0:N]
        c = conformal_reseries(b, 1.0, N)
        # Spot check against direct evaluation: pick a small w, compute both
        # B(t(w)) directly via 1/(1+t(w)), and the truncated polynomial in w.
        w0 = 0.1
        t0 = inverse_conformal_map(w0; a = 1.0)
        direct = 1 / (1 + t0)
        approx = sum(c[k+1] * w0^k for k in 0:N)
        @test isapprox(approx, direct; atol = 1e-6)
    end

    @testset "conformal_reseries Complex eltype" begin
        # Same setup but with complex `b`. Composition is linear in b, so
        # multiplying b by a complex constant scales the result identically.
        N = 12
        b_real = Float64[(-1.0)^k for k in 0:N]
        b_cplx = ComplexF64.(b_real) .* (1 + 0.1im)
        c = conformal_reseries(b_cplx, 1.0, N)
        @test eltype(c) === ComplexF64
        # Compare with the real reseries scaled by (1 + 0.1im).
        c_real = conformal_reseries(b_real, 1.0, N)
        @test c ≈ c_real .* (1 + 0.1im)
    end
end
