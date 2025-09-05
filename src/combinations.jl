"""
    combine_dempster(X, Y)

Combine two BPAs using Dempster's Rule of Combination.

See also: [`bpa`](@ref)
"""
function combine_dempster(X::BPA, Y::BPA)
    # if sum(values(X)) + sum(values(Y)) != one(Real) + one(Real)
    #     @warn "Mass assignments are not properly distributed!"
    # end

    # calculate the cross product of both mass assignments
    ps = Iterators.product(X, Y)

    # get all focal elements
    es = (keys(X) ∪ keys(Y))

    # subtract once so we don't have to do it for every iteration below
    one_minus_K = 1 - sum(p[1].second * p[2].second for p in ps if isempty(p[1].first ∩ p[2].first))

    # equalize BPA data types
    p_type = promote_type(first(eltype(X).types), first(eltype(Y).types))
    e_type = promote_type(last(eltype(X).types), last(eltype(Y).types))

    r = BPA{p_type, e_type}([
        e => sum(p[1].second * p[2].second for p in ps if (p[1].first ∩ p[2].first) == (e ∩ e); init = 0) / one_minus_K
        for e in es
    ])
    
    return r
end