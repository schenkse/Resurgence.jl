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

    @testset "pair map (singularities at ±i·a)" begin
        @testset "round-trip on real and complex t" begin
            for t in (0.1, 0.5, 1.0, 5.0)
                v = conformal_map_pair(t; a = 1)
                @test inverse_conformal_map_pair(v; a = 1) ≈ t
            end
            for t in (0.5 + 0.1im, 0.2 - 0.3im, 1.0 + 1.0im)
                v = conformal_map_pair(t; a = 1)
                @test inverse_conformal_map_pair(v; a = 1) ≈ t
            end
        end

        @testset "v(0) = 0 with real input stays real" begin
            @test conformal_map_pair(0.0; a = 1.0) == 0.0
            @test conformal_map_pair(0.0; a = 1.0) isa Real
            @test conformal_map_pair(0.5; a = 1.0) isa Real
        end

        @testset "argument validation" begin
            @test_throws ArgumentError conformal_map_pair(1.0; a = -1.0)
            @test_throws ArgumentError inverse_conformal_map_pair(0.5; a = 0.0)
            @test_throws ArgumentError conformal_reseries_pair([1.0], 1.0, -1)
        end

        @testset "singularities map to the unit circle (v = ±i)" begin
            # Map sends t = ±i·a to v = ±i on the boundary |v| = 1.
            v_plus = conformal_map_pair(im; a = 1.0)   # t = +i·a, a=1
            v_minus = conformal_map_pair(-im; a = 1.0) # t = -i·a
            @test isapprox(v_plus, im; atol = 1e-12)
            @test isapprox(v_minus, -im; atol = 1e-12)
            @test isapprox(abs(v_plus), 1; atol = 1e-12)
        end

        @testset "conformal_reseries_pair reproduces series at small v" begin
            # B(t) = 1/(1+t²) has singularities at t = ±i. Coefficients
            # are b_{2k} = (-1)^k, b_{2k+1} = 0.
            N = 12
            b = Float64[isodd(k) ? 0.0 : (-1.0)^(k ÷ 2) for k in 0:N]
            c = conformal_reseries_pair(b, 1.0, N)
            v0 = 0.2
            t0 = inverse_conformal_map_pair(v0; a = 1.0)
            direct = 1 / (1 + t0^2)
            approx = sum(c[k+1] * v0^k for k in 0:N)
            @test isapprox(approx, direct; atol = 1e-6)
        end

        @testset "Complex eltype propagates through conformal_reseries_pair" begin
            N = 12
            b_real = Float64[isodd(k) ? 0.0 : (-1.0)^(k ÷ 2) for k in 0:N]
            b_cplx = ComplexF64.(b_real) .* (1 + 0.1im)
            c = conformal_reseries_pair(b_cplx, 1.0, N)
            @test eltype(c) === ComplexF64
            c_real = conformal_reseries_pair(b_real, 1.0, N)
            @test c ≈ c_real .* (1 + 0.1im)
        end
    end
end
