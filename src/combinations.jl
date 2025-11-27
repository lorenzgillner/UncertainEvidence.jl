"""
    combine_dempster(X, Y)

Combine two BPAs using Dempster's Rule of Combination.

See also: [`combine_yager`](@ref), [`bpa`](@ref).
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

    r = BPA{p_type,e_type}(
        (
            e => sum(p[1].second * p[2].second for p in ps if (p[1].first ∩ p[2].first) == (e ∩ e); init=0) / one_minus_K
            for e in es
        )...
    )

    return r
end

"""
    combine_yager(X::BPA, Y::BPA)

Combine two BPAs using Yager's Rule of Combination.

See also: [`combine_dempster`](@ref), [`bpa`](@ref).
"""
function combine_yager(X::BPA, Y::BPA)
    # calculate the cross product of both mass assignments
    ps = collect(Iterators.product(collect(X), collect(Y)))

    # get all focal elements
    es = (keys(X) ∪ keys(Y))
    
    # compute K, the mass of conflict
    K = sum(p[1].second * p[2].second for p in ps if isempty(p[1].first ∩ p[2].first))
    
    r = BPA([
        e => sum(p[1].second * p[2].second for p in ps if (p[1].first ∩ p[2].first) == (e ∩ e); init = 0)
        for e in es
    ])

    all_keys = reduce(union, keys(r))
    
    r[all_keys] = r[all_keys] + K
    
    return r
end