using Test
using Resurgence: obtain_poles_real, move_poles, poles_regularized
using Polynomials: Polynomial, coeffs

@testset "poles" begin
    @testset "obtain_poles_real" begin
        # (t - 2)(t + 3)(t² + 1) = t⁴ + t³ - 5t² + t - 6  (constant first: -6, 1, -5, 1, 1)
        p = Polynomial([-6.0, 1.0, -5.0, 1.0, 1.0])
        rr = obtain_poles_real(coeffs(p))
        @test length(rr) == 2
        @test 2.0 in round.(rr; digits = 6)
        @test -3.0 in round.(rr; digits = 6)
    end

    @testset "move_poles shifts only positive real roots" begin
        roots = ComplexF64[2.0, -1.0, 0.5 + 0.3im]
        out = move_poles(roots, 1e-3)
        @test eltype(out) === ComplexF64
        @test imag(out[1]) ≈ 1e-3      # shifted
        @test out[2] == ComplexF64(-1.0, 0.0)        # not shifted
        @test out[3] == ComplexF64(0.5, 0.3)         # not shifted
    end

    @testset "poles_regularized normalizes constant term" begin
        # q(t) = 1 - t  has root t = 1 on the positive real axis.
        qr = poles_regularized([1.0, -1.0], 1e-6)
        c = coeffs(qr)
        @test c[1] ≈ 1.0
        # Denominator should now have a complex root with small +iε imaginary part.
        @test abs(qr(1.0)) > 0   # no longer vanishes at t=1
    end
end
