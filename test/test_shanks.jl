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
end
