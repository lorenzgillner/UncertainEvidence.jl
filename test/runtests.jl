using UncertainEvidence
using Test
using LazySets
using LinearAlgebra

# ==================================
# Three Color Example from Wikipedia
# ==================================
# See: https://en.wikipedia.org/wiki/Dempster%E2%80%93Shafer_theory#Bayesian_approximation

# First sensor
m1 = Dict(
    Set("r") => 0.35,
    Set("y") => 0.25,
    Set("g") => 0.15,
    Set("ry") => 0.06,
    Set("rg") => 0.05,
    Set("yg") => 0.04,
    Set("ryg") => 0.1
)

# Second sensor; notice how this one is missing the mass for Ω:
m2 = Dict(
    Set("r") => 0.11,
    Set("y") => 0.21,
    Set("g") => 0.33,
    Set("ry") => 0.21,
    Set("rg") => 0.01,
    Set("yg") => 0.03,
)

# Combined data, rounded to two decimal places
m12 = Dict(
    Set(['r'])           => 0.32,
    Set(['g', 'r'])      => 0.01,
    Set(['y'])           => 0.33,
    Set(['g', 'y'])      => 0.01,
    Set(['g', 'y', 'r']) => 0.02,
    Set(['g'])           => 0.24,
    Set(['y', 'r'])      => 0.07
)

mc = combine_dempster(m1, m2)

@test sum(values(dss(m2))) == 1.0
@test keys(mc) == keys(m12)
@test round.(values(mc), digits=2) == collect(values(m12))

# ==================
# Earthquake Example
# ==================
# Inspired by Wang & Klir: "Fuzzy Measure Theory"

# Epicenter of the earthquake
B = Ball2([2.0, 1.0], 1.0)

# Estimates for the earthquake's epicenter
E1 = Ball2([2.5, 0.75], 0.25)
E2 = Ball2([1.8, 1.8], 0.5)
E3 = Ball2([2.5, 2.5], 0.25)
E4 = Ball2([2.7, 2.5], 0.2)

estimates = [E1, E2, E3, E4]
masses = fill(1/4, 4)
me = Dict(zip(estimates, masses))

@test bel(B, me) == 0.25
@test pls(B, me) == 0.5