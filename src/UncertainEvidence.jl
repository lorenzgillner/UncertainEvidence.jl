module UncertainEvidence

import Base: @warn,
    ∩, ∪, ∈, ∉, +, -, *, /,
    Dict, Pair, Tuple, zero, one,
    keys, values, eltype,
    reduce, sum, first, last,
    in, intersect,
    issubset, isdisjoint, isempty
    Iterators.product

import Combinatorics: combinations

export
    BPA, bpa,
    focalelements, masses, omega,
    redistribute!,
    bel, pls
include("dst.jl")

export combine_dempster, combine_yager
include("combinations.jl")

end
