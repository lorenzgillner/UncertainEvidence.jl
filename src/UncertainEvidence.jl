module UncertainEvidence

include("dst.jl")
include("combinations.jl")

export
    BPA, bpa,
    redistribute!,
    bel, pls

export combine_dempster

# TODO p-boxes

end
