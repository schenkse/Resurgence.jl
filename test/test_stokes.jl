using Test
using Resurgence
using SpecialFunctions: gamma

# Build coefficients with the exact asymptotic
#
#     a[k+1] = A · Γ(k+β) / S^(k+β) · (1 + c₁/k + c₂/k² + …)
#
# for k = 0, 1, …, N-1. For k = 0 the 1/k correction is undefined; the test
# series themselves use only the leading Γ(k+β)/S^(k+β) form (no c terms).
function _exact_asym(::Type{T}, N, S, β, A) where {T<:Number}
    out = Vector{T}(undef, N)
    @inbounds for k in 0:N-1
        out[k+1] = T(A) * gamma(T(k) + T(β)) / T(S)^(T(k) + T(β))
    end
    return out
end

@testset "stokes — large-order extraction" begin
    @testset "unsigned shifted Stieltjes a_k = (k+1)! ⇒ S=1, β=2, A=1" begin
        # a[k+1] = (k+1)!  =  Γ(k+2),  matches S=1, β=2, A=1.
        a = BigFloat[BigFloat(factorial(big(k + 1))) for k in 0:49]
        @test isapprox(stokes_action(a),   BigFloat(1); atol = BigFloat("1e-15"))
        @test isapprox(stokes_exponent(a, BigFloat(1)), BigFloat(2); atol = BigFloat("1e-15"))
        @test isapprox(stokes_constant(a, BigFloat(1), BigFloat(2)),
                       BigFloat(1); atol = BigFloat("1e-15"))

        fit = stokes_fit(a)
        @test isapprox(fit.S, BigFloat(1); atol = BigFloat("1e-15"))
        @test isapprox(fit.β, BigFloat(2); atol = BigFloat("1e-15"))
        @test isapprox(fit.A, BigFloat(1); atol = BigFloat("1e-15"))
        @test isempty(fit.c)
    end

    @testset "alternating shifted Stieltjes a_k = (-1)^k (k+1)! ⇒ S=-1, β=2, A=1" begin
        a = BigFloat[(-1)^k * BigFloat(factorial(big(k + 1))) for k in 0:49]
        @test isapprox(stokes_action(a),   BigFloat(-1); atol = BigFloat("1e-15"))
        @test isapprox(stokes_exponent(a, BigFloat(-1)), BigFloat(2); atol = BigFloat("1e-15"))
        @test isapprox(stokes_constant(a, BigFloat(-1), BigFloat(2)),
                       BigFloat(1); atol = BigFloat("1e-15"))
    end

    @testset "plain Stieltjes a_k = (-1)^k k! ⇒ S=-1, β=1, A=-1" begin
        # a[k+1] = (-1)^k Γ(k+1).  Match A·Γ(k+β)/S^(k+β):  β = 1, S = -1,
        # and A = (-1)^k · (-1)^{k+1} = -1 (the leftover (-1)^{2k+1}).
        a = BigFloat[(-1)^k * BigFloat(factorial(big(k))) for k in 0:49]
        @test isapprox(stokes_action(a), BigFloat(-1); atol = BigFloat("1e-12"))
        @test isapprox(stokes_exponent(a, BigFloat(-1)), BigFloat(1); atol = BigFloat("1e-12"))
        @test isapprox(stokes_constant(a, BigFloat(-1), BigFloat(1)),
                       BigFloat(-1); atol = BigFloat("1e-12"))
    end

    @testset "subleading c-coefficient recovery" begin
        # Construct a[k+1] = Γ(k+β)/S^(k+β) · (1 + c₁/k + c₂/k²) for k ≥ 1
        # and a[1] = limit (anything finite — only the tail matters).
        S, β, A = BigFloat(1), BigFloat(2), BigFloat(1)
        c1, c2 = BigFloat("0.5"), BigFloat("-0.25")
        N = 60
        a = _exact_asym(BigFloat, N, S, β, A)
        @inbounds for k in 1:N-1
            a[k+1] *= (BigFloat(1) + c1 / BigFloat(k) + c2 / BigFloat(k)^2)
        end
        # Use larger depth than default to reach asymptotic regime.
        fit = stokes_fit(a; subleading = 2)
        @test isapprox(fit.S, S; atol = BigFloat("1e-6"))
        @test isapprox(fit.β, β; atol = BigFloat("1e-6"))
        @test isapprox(fit.A, A; atol = BigFloat("1e-6"))
        @test length(fit.c) == 2
        @test isapprox(fit.c[1], c1; atol = BigFloat("1e-3"))
        @test isapprox(fit.c[2], c2; atol = BigFloat("1e-2"))
    end

    @testset "Float64 element-type genericity" begin
        a = Float64[Float64(factorial(big(k + 1))) for k in 0:24]
        @test isapprox(stokes_action(a), 1.0; atol = 1e-6)
        @test isapprox(stokes_exponent(a, 1.0), 2.0; atol = 1e-6)
        @test isapprox(stokes_constant(a, 1.0, 2.0), 1.0; atol = 1e-6)
        fit = stokes_fit(a)
        @test fit.S isa Float64
        @test fit.β isa Float64
        @test fit.A isa Float64
    end

    @testset "Complex coefficients propagate complex S, β, A" begin
        # alternating shifted: real S=-1, real β=2, real A=1, but as complex.
        a = ComplexF64[ComplexF64((-1)^k * factorial(big(k + 1))) for k in 0:24]
        S = stokes_action(a)
        @test S isa Complex
        @test isapprox(S, ComplexF64(-1); atol = 1e-6)
    end

    @testset "argument validation" begin
        @test_throws ArgumentError stokes_action(Float64[1.0, 2.0])  # length < 4
        @test_throws ArgumentError stokes_fit(Float64[1.0, 1.0, 1.0, 1.0]; subleading = -1)
        # zero in the middle of the ratio sequence
        a_bad = Float64[1.0, 0.0, 1.0, 1.0, 1.0]
        @test_throws ArgumentError stokes_action(a_bad)
    end
end
