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
    BPA, bpa, BaseType,
    focalelements, masses, omega,
    redistribute!
include("bpa.jl")

export bel, pls
include("dst.jl")

export combine_dempster, combine_yager
include("combinations.jl")

end
