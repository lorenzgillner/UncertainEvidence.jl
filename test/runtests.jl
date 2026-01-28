using UncertainEvidence
using Test

@testset verbose = true "UncertainEvidence" begin
    @testset "BPA creation" begin
        @testset "From Dict" begin
            @testset "With explicit Ω" begin
                d = Dict(Set([:a]) => 0.5, Set([:b]) => 0.5, Set([:a, :b]) => 0.0)
                Ω = Set([:a, :b])
                X = BPA(d, Ω)

                @test X[:a] == 0.5
                @test X[:b] == 0.5
                @test X[:a, :b] == 0.0
                @test length(X) == 3
                @test focalelements(X) == Set([Set([:a]), Set([:b]), Set([:a, :b])])
                @test frame(X) == Set([:a, :b])
            end

            @testset "With explicit Ω (as kwarg)" begin
                d = Dict(Set([:a]) => 0.5, Set([:b]) => 0.5)
                X = BPA(d, Ω=Set([:a, :b]))

                @test X[:a] == 0.5
                @test X[:b] == 0.5
                @test X[:a, :b] == 0.0
                @test length(X) == 3
                @test focalelements(X) == Set([Set([:a]), Set([:b]), Set([:a, :b])])
                @test frame(X) == Set([:a, :b])
            end

            @testset "With implicit Ω" begin
                d = Dict(Set([:a]) => 0.5, Set([:b]) => 0.5)
                X = BPA(d)

                @test X[:a] == 0.5
                @test X[:b] == 0.5
                @test X[:a, :b] == 0.0
                @test length(X) == 3
                @test focalelements(X) == Set([Set([:a]), Set([:b]), Set([:a, :b])])
                @test frame(X) == Set([:a, :b])
            end

            @testset "Splat operator" begin
                d = Dict(Set([:a]) => 0.5, Set([:b]) => 0.5)
                X = BPA(d...)

                @test X[:a] == 0.5
                @test X[:b] == 0.5
                @test X[:a, :b] == 0.0
                @test length(X) == 3
                @test focalelements(X) == Set([Set([:a]), Set([:b]), Set([:a, :b])])
                @test frame(X) == Set([:a, :b])
            end

            @testset "Non-set types" begin
                d = Dict(:a => 0.5, :b => 0.5)
                X = BPA(d)

                @test X[:a] == 0.5
                @test X[:b] == 0.5
                @test X[:a, :b] == 0.0
                @test length(X) == 3
                @test focalelements(X) == Set([Set([:a]), Set([:b]), Set([:a, :b])])
                @test frame(X) == Set([:a, :b])
            end
        end

        @testset "From Pairs" begin
            @testset "With explicit Ω" begin
                X = BPA(:a => 0.5, :b => 0.5; Ω=Set([:a, :b]))

                @test X[:a] == 0.5
                @test X[:b] == 0.5
                @test X[:a, :b] == 0.0
                @test length(X) == 3
                @test focalelements(X) == Set([Set([:a]), Set([:b]), Set([:a, :b])])
                @test frame(X) == Set([:a, :b])
            end

            @testset "With implicit Ω" begin
                X = BPA(:a => 0.5, :b => 0.5)

                @test X[:a] == 0.5
                @test X[:b] == 0.5
                @test X[:a, :b] == 0.0
                @test length(X) == 3
                @test focalelements(X) == Set([Set([:a]), Set([:b]), Set([:a, :b])])
                @test frame(X) == Set([:a, :b])
            end
        end
    end

    @testset "Combination rules" begin
        @testset "Dempster's rule" begin
            X₁ = BPA(:a => 0.5, :b => 0.5)
            @test frame(X₁) == Set([:a, :b])
            @test X₁[:a, :b] == 0.0

            X₂ = BPA(:a => 0.25, :b => 0.75)
            @test frame(X₂) == Set([:a, :b])
            @test X₂[:a, :b] == 0.0

            X = combine_dempster(X₁, X₂)

            @test isnormal(X)

            @test X[:a] == 0.25
            @test X[:b] == 0.75
            @test X[:a, :b] == 0.0
        end

        # @testset "Yager's rule" begin
        # end
    end

    @testset "Focal element types" begin
        @testset "Characters" begin
            # test combination rules, based on Zadeh's paradox, see:
            # https://doi.org/10.1609/aimag.v5i3.452

            X₁ = BPA(
                'A' => 0.99,
                'B' => 0.01,
                'C' => 0.00,
            )

            @test isnormal(X₁)

            X₂ = BPA(
                'A' => 0.00,
                'B' => 0.01,
                'C' => 0.99,
            )

            @test isnormal(X₂)

            X = combine_dempster(X₁, X₂)

            @test X['A'] == 0.0
            @test X['B'] ≈ 1.0
            @test X['C'] == 0.0
        end

        @testset "Sets of characters" begin
            # three colors example from Wikipedia, see:
            # https://en.wikipedia.org/wiki/Dempster%E2%80%93Shafer_theory#Bayesian_approximation

            # first sensor
            m₁ = BPA(
                Set("r") => 0.35,
                Set("y") => 0.25,
                Set("g") => 0.15,
                Set("ry") => 0.06,
                Set("rg") => 0.05,
                Set("yg") => 0.04,
                Set("ryg") => 0.1,
            )

            @test isnormal(m₁)

            # second sensor; notice how this one is missing the mass assignment for Ω:
            m₂ = BPA(
                Set("r") => 0.11,
                Set("y") => 0.21,
                Set("g") => 0.33,
                Set("ry") => 0.21,
                Set("rg") => 0.01,
                Set("yg") => 0.03,
            )

            @test isnormal(m₂)

            # combined data, rounded to two decimal places
            m₁₂ = BPA(
                Set("r") => 0.32,
                Set("y") => 0.33,
                Set("g") => 0.24,
                Set("yr") => 0.07,
                Set("gr") => 0.01,
                Set("gy") => 0.01,
                Set("gyr") => 0.02,
            )

            @test isnormal(m₁₂)

            mc = combine_dempster(m₁, m₂)

            @test isnormal(mc)

            @test focalelements(mc) == focalelements(m₁₂)

            @test sort(round.(masses(mc), digits=2)) == sort(collect(masses(m₁₂)))
        end

        @testset "Sets of strings" begin
            # Zadeh's paradox, but this time with more descriptive focal elements

            X₁ = BPA(
                "concussion" => 0.99,
                "tumor" => 0.01,
                "migraine" => 0.00,
            )

            @test isnormal(X₁)

            X₂ = BPA(
                "concussion" => 0.00,
                "tumor" => 0.01,
                "migraine" => 0.99,
            )

            @test isnormal(X₂)

            X = combine_dempster(X₁, X₂)

            @test isnormal(X)

            @test X["concussion"] == 0.0
            @test X["tumor"] ≈ 1.0
            @test X["migraine"] == 0.0
        end
    end

    # @testset "Ellipsoids (ℝ²)" begin
    # 	# Earthquake example, inspired by:
    # 	# Z. Wang, G. J. Klir (2013): "Fuzzy measure theory"

    # 	# Epicenter of the earthquake
    # 	B = Ball2([2.0, 1.0], 1.0)

    # 	# Estimates for the earthquake's epicenter
    # 	E1 = Ball2([2.5, 0.75], 0.25)
    # 	E2 = Ball2([1.8, 1.8], 0.5)
    # 	E3 = Ball2([2.5, 2.5], 0.25)
    # 	E4 = Ball2([2.7, 2.5], 0.2)

    # 	estimates = [E1, E2, E3, E4]
    # 	masses = fill(1.0 / 4, 4)
    # 	me = BPA(zip(estimates, masses))

    # 	@test bel(B, me) == 0.25
    # 	@test pls(B, me) == 0.5
    # end
end
