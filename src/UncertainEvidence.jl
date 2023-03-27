module UncertainEvidence

include("dst.jl")
include("combinations.jl")

export
    DSS, dss,
    redistribute!,
    bel, pls

export combine_dempster

# TODO p-boxes

end
