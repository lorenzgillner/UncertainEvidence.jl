using UncertainEvidence

using Test
using LazySets
using LinearAlgebra

@testset verbose = true "UncertainEvidence" begin
    @testset "basics" begin
        # test correct redistribution
        A = bpa(
            Set("a") => 0.1,
            Set("b") => 0.2
        )
        @test A[Set("ab")] == 0.7
    end

    @testset verbose = true "focal element types" begin
        @testset "characters" begin
            # test combination rules, based on Zadeh's paradox, see:
            # https://doi.org/10.1609/aimag.v5i3.452

            X1 = BPA(
                'A' => 0.99,
                'B' => 0.01,
                'C' => 0.00
            )

            X2 = BPA(
                'A' => 0.00,
                'B' => 0.01,
                'C' => 0.99
            )

            X12 = combine_dempster(X1, X2)

            @test X12['A'] == 0.0
            @test X12['B'] ≈ 1.0
            @test X12['C'] == 0.0
        end

        @testset "sets of characters" begin
            # three colors example from Wikipedia, see:
            # https://en.wikipedia.org/wiki/Dempster%E2%80%93Shafer_theory#Bayesian_approximation

            # first sensor
            m1 = BPA(
                Set("r") => 0.35,
                Set("y") => 0.25,
                Set("g") => 0.15,
                Set("ry") => 0.06,
                Set("rg") => 0.05,
                Set("yg") => 0.04,
                Set("ryg") => 0.1,
            )

            # second sensor; notice how this one is missing the mass assignment for Ω:
            m2 = BPA(
                Set("r") => 0.11,
                Set("y") => 0.21,
                Set("g") => 0.33,
                Set("ry") => 0.21,
                Set("rg") => 0.01,
                Set("yg") => 0.03,
            )

            # combined data, rounded to two decimal places
            m12 = BPA(
                Set("r") => 0.32,
                Set("y") => 0.33,
                Set("g") => 0.24,
                Set("yr") => 0.07,
                Set("gr") => 0.01,
                Set("gy") => 0.01,
                Set("gyr") => 0.02,
            )

            @test sum(values(redistribute!(m2))) == one(Real)

            mc = combine_dempster(m1, m2)

            @test keys(mc) == keys(m12)

            @test sort(round.(values(mc), digits=2)) == sort(collect(values(m12)))
        end

        @testset "sets of strings" begin
            # Zadeh's paradox again, but this time with more descriptive focal elements

            X1 = BPA(
                ["concussion"] => 0.99,
                ["tumor"] => 0.01,
                ["migraine"] => 0.00
            )

            X2 = BPA(
                ["concussion"] => 0.00,
                ["tumor"] => 0.01,
                ["migraine"] => 0.99
            )

            X12 = combine_dempster(X1, X2)

            @test X12[["concussion"]] == 0.0
            @test X12[["tumor"]] ≈ 1.0
            @test X12[["migraine"]] == 0.0
        end

        @testset "balls (ℝ²)" begin
            # earthquake example, inspired by:
            # Z. Wang, G. J. Klir (2013): "Fuzzy measure theory"

            # epicenter of the earthquake
            B = Ball2([2.0, 1.0], 1.0)

            # estimates for the earthquake's epicenter
            E1 = Ball2([2.5, 0.75], 0.25)
            E2 = Ball2([1.8, 1.8], 0.5)
            E3 = Ball2([2.5, 2.5], 0.25)
            E4 = Ball2([2.7, 2.5], 0.2)

            estimates = [E1, E2, E3, E4]
            masses = fill(1.0 / 4, 4)
            me = bpa(zip(estimates, masses))

            @test bel(B, me) == 0.25
            @test pls(B, me) == 0.5
        end
    end
end