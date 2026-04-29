using Test
using Resurgence: inv_factorials, chop!, sparsify!, split_vector

@testset "utils" begin
    @testset "inv_factorials Float64" begin
        f = inv_factorials(Float64, 6)
        @test f == [1.0, 1.0, 1/2, 1/6, 1/24, 1/120]
    end

    @testset "inv_factorials BigFloat past Float64 overflow" begin
        # length 30 is well past where factorial(::Int) overflows (k = 21).
        f = inv_factorials(BigFloat, 30)
        @test f[1] == big"1"
        @test f[2] == big"1"
        # spot check k = 25: 1/25! computed independently. Iterative division
        # accumulates a few ulps of BigFloat rounding — well within precision.
        @test isapprox(f[26], 1 / factorial(big(25)); rtol = BigFloat("1e-70"))
        @test eltype(f) === BigFloat
    end

    @testset "inv_factorials zero length" begin
        @test isempty(inv_factorials(Float64, 0))
    end

    @testset "chop! real" begin
        x = [1e-20, 1.0, -2e-20, 3.0]
        chop!(x, 1e-15)
        @test x == [0.0, 1.0, 0.0, 3.0]
        @test eltype(x) === Float64
    end

    @testset "chop! complex preserves eltype" begin
        x = ComplexF64[1e-20 + 1.0im, 1.0 + 1e-20im, 1e-20 + 1e-20im]
        chop!(x, 1e-15)
        @test x == ComplexF64[0.0 + 1.0im, 1.0 + 0.0im, 0.0 + 0.0im]
        @test eltype(x) === ComplexF64
    end

    @testset "sparsify! aliases chop!" begin
        x = [1e-20, 1.0]
        sparsify!(x, 1e-15)
        @test x == [0.0, 1.0]
    end

    @testset "split_vector by Int" begin
        v = [1, 2, 3, 4, 5]
        out = split_vector(v, 2)
        @test out == [[1, 2], [3, 4, 5]]
    end

    @testset "split_vector by lengths" begin
        v = collect(1:6)
        out = split_vector(v, [1, 2, 3])
        @test out == [[1], [2, 3], [4, 5, 6]]
    end
end
