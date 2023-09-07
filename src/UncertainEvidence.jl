__precompile__(true)

module UncertainEvidence

import Base: @warn,
    ∩, ∪, ∈, ∉, +, -, *, /,
    Dict, Pair, Tuple, zero, one,
    keys, values, eltype,
    reduce, sum, first, last,
    in, intersect,
    issubset, isdisjoint, isempty
    Iterators.product

export
    BPA, bpa,
    redistribute!,
    bel, pls
include("dst.jl")

export combine_dempster
include("combinations.jl")

end
