using Test
using Resurgence

@testset "TransSeries" begin
    @testset "Sector construction promotes mixed inputs" begin
        s = Sector(1, 0.0, [1.0, 2.0])
        @test s isa Sector{Float64}
        @test s.S === 1.0
        @test s.β === 0.0
        @test s.a == [1.0, 2.0]

        sc = Sector(1, 0.0, ComplexF64[1, 2 + im])
        @test sc isa Sector{ComplexF64}
        @test sc.S === 1.0 + 0.0im

        sb = Sector(big"1", big"0.5", [big"1.0", big"2.0"])
        @test sb isa Sector{BigFloat}
    end

    @testset "TransSeries(a::Vector) wraps as perturbative sector" begin
        a = [1.0, 2.0, 3.0]
        ts = TransSeries(a)
        @test ts isa TransSeries{Float64}
        @test length(ts.sectors) == 1
        @test ts.sectors[1].S === 0.0
        @test ts.sectors[1].β === 0.0
        @test ts.sectors[1].a == a
    end

    @testset "TransSeries from sector list" begin
        ts = TransSeries([Sector(0.0, 0.0, [1.0, 2.0]),
                          Sector(1.0, 0.5, [3.0])])
        @test ts isa TransSeries{Float64}
        @test length(ts.sectors) == 2

        # Empty list is rejected (no T to infer).
        @test_throws ArgumentError TransSeries(Sector[])
    end

    @testset "+ concatenates and merges equal (S, β)" begin
        ts1 = TransSeries([Sector(0.0, 0.0, [1.0]),
                           Sector(1.0, 0.0, [2.0])])
        ts2 = TransSeries([Sector(0.0, 0.0, [3.0, 4.0]),
                           Sector(2.0, 0.0, [5.0])])
        ts = ts1 + ts2
        @test length(ts.sectors) == 3
        # (0, 0): [1] + [3, 4] = [4, 4]
        idx0 = findfirst(s -> s.S == 0.0 && s.β == 0.0, ts.sectors)
        @test ts.sectors[idx0].a == [4.0, 4.0]
        # (1, 0): unchanged
        idx1 = findfirst(s -> s.S == 1.0 && s.β == 0.0, ts.sectors)
        @test ts.sectors[idx1].a == [2.0]
        # (2, 0): unchanged
        idx2 = findfirst(s -> s.S == 2.0 && s.β == 0.0, ts.sectors)
        @test ts.sectors[idx2].a == [5.0]
    end

    @testset "+ keeps sectors with different β separate" begin
        ts1 = TransSeries([Sector(0.0, 0.0, [1.0])])
        ts2 = TransSeries([Sector(0.0, 1.0, [2.0])])
        ts = ts1 + ts2
        @test length(ts.sectors) == 2
    end

    @testset "Unary - and binary - work" begin
        ts = TransSeries([Sector(0.0, 0.0, [1.0, 2.0]),
                          Sector(1.0, 0.0, [3.0])])
        @test (-ts).sectors[1].a == [-1.0, -2.0]
        @test (-ts).sectors[2].a == [-3.0]
        @test (ts - ts).sectors[1].a == [0.0, 0.0]
    end

    @testset "scalar * scales each sector" begin
        ts = TransSeries([Sector(0.0, 0.0, [1.0, 2.0]),
                          Sector(1.0, 0.0, [3.0])])
        ts2 = 2 * ts
        @test ts2.sectors[1].a == [2.0, 4.0]
        @test ts2.sectors[2].a == [6.0]
        # Commutative form is provided.
        @test (ts * 0.5).sectors[1].a == [0.5, 1.0]
    end

    @testset "TransSeries * TransSeries: single-sector Cauchy product" begin
        # (1, 0, [1, 2]) * (1, 0, [1, 3]) → (2, 0, [1, 5, 6])
        s1 = TransSeries([Sector(1.0, 0.0, [1.0, 2.0])])
        s2 = TransSeries([Sector(1.0, 0.0, [1.0, 3.0])])
        prod = s1 * s2
        @test length(prod.sectors) == 1
        @test prod.sectors[1].S === 2.0
        @test prod.sectors[1].β === 0.0
        @test prod.sectors[1].a == [1.0, 5.0, 6.0]
    end

    @testset "TransSeries * TransSeries: cross-products merge by (S, β)" begin
        # ts = sector_A + sector_B; ts * ts produces 4 raw products but
        # A*B and B*A share (S, β) and must merge.
        sA = Sector(1.0, 0.0, [1.0])  # e^{-1/g}·1
        sB = Sector(2.0, 0.0, [1.0])  # e^{-2/g}·1
        ts = TransSeries([sA, sB])
        ts2 = ts * ts
        # Expected sectors: (2, 0, [1]), (3, 0, [1+1]=[2]), (4, 0, [1])
        @test length(ts2.sectors) == 3
        actions = sort([s.S for s in ts2.sectors])
        @test actions == [2.0, 3.0, 4.0]
        idx3 = findfirst(s -> s.S == 3.0, ts2.sectors)
        @test ts2.sectors[idx3].a == [2.0]
    end

    @testset "transseries_exp at small orders" begin
        ts = TransSeries([Sector(1.0, 0.0, [1.0])])

        # order = 0 ⇒ just 1
        ex0 = transseries_exp(ts; order = 0)
        @test length(ex0.sectors) == 1
        @test ex0.sectors[1].S === 0.0
        @test ex0.sectors[1].β === 0.0
        @test ex0.sectors[1].a == [1.0]

        # order = 1 ⇒ 1 + ts
        ex1 = transseries_exp(ts; order = 1)
        @test length(ex1.sectors) == 2
        idx0 = findfirst(s -> s.S == 0.0, ex1.sectors)
        idx1 = findfirst(s -> s.S == 1.0, ex1.sectors)
        @test ex1.sectors[idx0].a == [1.0]
        @test ex1.sectors[idx1].a == [1.0]

        # order = 3 ⇒ dilute instanton gas truncated at 3 instantons
        # ts^k = (k·S, 0, [1])  so coefficient of e^{-k·S/g} is 1/k!
        ex3 = transseries_exp(ts; order = 3)
        actions = sort([s.S for s in ex3.sectors])
        @test actions == [0.0, 1.0, 2.0, 3.0]
        for (k, want) in [(0, 1.0), (1, 1.0), (2, 0.5), (3, 1/6)]
            idx = findfirst(s -> s.S == Float64(k), ex3.sectors)
            @test ex3.sectors[idx].a == [want]
        end
    end

    @testset "transseries_exp rejects negative order" begin
        ts = TransSeries([Sector(1.0, 0.0, [1.0])])
        @test_throws ArgumentError transseries_exp(ts; order = -1)
    end

    @testset "BigFloat eltype round-trips" begin
        ts = TransSeries([Sector(big"1.0", big"0.0", BigFloat[1, 2])])
        @test ts isa TransSeries{BigFloat}
        @test (ts + ts) isa TransSeries{BigFloat}
        @test (ts * ts) isa TransSeries{BigFloat}
        ex = transseries_exp(ts; order = 2)
        @test ex isa TransSeries{BigFloat}
        # 1 + ts + ts²/2 has sector at action 2 with value [1, 4, 4]/2 = [0.5, 2, 2]
        idx2 = findfirst(s -> s.S == big"2.0", ex.sectors)
        @test ex.sectors[idx2].a ≈ BigFloat[0.5, 2.0, 2.0]
    end

    @testset "ComplexF64 eltype round-trips" begin
        ts = TransSeries([Sector(1.0 + 0.0im, 0.0 + 0.0im, ComplexF64[1, im])])
        @test ts isa TransSeries{ComplexF64}
        prod = ts * ts
        @test prod isa TransSeries{ComplexF64}
        # [1, i] ⋆ [1, i] = [1, 2i, -1]
        @test prod.sectors[1].a ≈ ComplexF64[1, 2im, -1]
    end

    @testset "Mixed eltype + promotes" begin
        ts_f = TransSeries(Float64[1.0, 2.0])
        ts_b = TransSeries(BigFloat[3.0, 4.0])
        ts = ts_f + ts_b
        @test ts isa TransSeries{BigFloat}
        @test ts.sectors[1].a ≈ BigFloat[4.0, 6.0]
    end
end
