module UncertainEvidence

import Base: @warn,
    ∩, ∪, ∈, ∉, +, -, *, /,
    Dict, Pair, Tuple, zero, one,
    keys, values, eltype,
    reduce, sum, first, last,
    in, intersect,
    issubset, isdisjoint, isempty
    Iterators.product

include("bpa.jl")
export
    BPA, BaseType,
    focalelements, masses, totalmass, frame

include("dst.jl")
export bel, pls, isnormal, normalize!

include("combinations.jl")
export combine_dempster, combine_yager

end # module