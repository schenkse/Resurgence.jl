using Test
using Resurgence

@testset "borel transforms" begin
    @testset "borel_transform basic" begin
        a = Float64[1, 2, 4, 6, 8]
        B = borel_transform(a)
        @test B ≈ [1.0, 2.0, 2.0, 1.0, 1/3]
    end

    @testset "borel_transform handles length > 21 in Float64" begin
        # Stieltjes coefficients (-1)^k k! truncate Float64 at large k, but the
        # transform B[k] = (-1)^k stays exact for all k.
        N = 30
        a = Float64[(-1.0)^k * factorial(big(k)) for k in 0:N-1]
        # convert big factorial back to Float64 then transform
        af = Float64.(a)
        # note: factorial(big(20)) > 1.0e18 so af loses precision past k=20.
        # We only test the leading entries.
        B = borel_transform(af)
        @test B[1] ≈ 1.0
        @test B[2] ≈ -1.0
        @test B[3] ≈ 1.0
        @test B[4] ≈ -1.0
    end

    @testset "borel_transform BigFloat preserves type" begin
        a = BigFloat[(-1)^k * factorial(big(k)) for k in 0:25]
        B = borel_transform(a)
        @test eltype(B) === BigFloat
        @test all(B[1:25] .≈ BigFloat[(-1)^k for k in 0:24])
    end

    @testset "borel_leroy_transform reduces to borel at b=0" begin
        a = [1.0, 2.0, 4.0, 6.0, 8.0]
        Bb = borel_leroy_transform(a, 0.0)
        B  = borel_transform(a)
        @test Bb ≈ B
    end

    @testset "borel_ratios" begin
        b = [1.0, 2.0, 6.0, 24.0]
        r = borel_ratios(b)
        @test r ≈ [2.0, 3.0, 4.0]
    end

    @testset "borel_ratios preserves BigFloat" begin
        b = BigFloat[1, 2, 6, 24]
        r = borel_ratios(b)
        @test eltype(r) === BigFloat
    end

    @testset "borel_ratios rejects zero divisors" begin
        @test_throws ArgumentError borel_ratios([1.0, 0.0, 1.0])
        @test_throws ArgumentError borel_ratios([0.0, 1.0])
    end

    @testset "Complex eltype" begin
        # Type-genericity smoke test: complex `a` should round-trip through
        # the Borel transforms preserving ComplexF64 and reducing correctly.
        a = ComplexF64[1 + 0im, 2 + 1im, 4 - 1im, 6 + 2im, 8]
        B = borel_transform(a)
        @test eltype(B) === ComplexF64
        # B[k] = a[k] / k!
        @test B ≈ a ./ Float64[factorial(k) for k in 0:length(a)-1]

        # borel_leroy_transform at b=0 must agree with borel_transform.
        Bb = borel_leroy_transform(a, 0.0)
        @test eltype(Bb) === ComplexF64
        @test Bb ≈ B

        # borel_ratios on complex input preserves ComplexF64.
        r = borel_ratios(B)
        @test eltype(r) === ComplexF64
    end
end
