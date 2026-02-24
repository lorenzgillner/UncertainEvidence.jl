"""
    combine_dempster(X, Y)

Combine two BPAs using Dempster's Rule of Combination.

See also: [`combine_yager`](@ref), [`bpa`](@ref).
"""
function combine_dempster(X::BPA{K,V}, Y::BPA{K,V}) where {K,V}
    # Calculate the cross product of both mass assignments
    ps = Iterators.product(X, Y)

    # Get all focal elements
    focal_elements = (focalelements(X) ∪ focalelements(Y))

    k = sum(p[1].second * p[2].second for p in ps if isempty(p[1].first ∩ p[2].first))

    # Subtract once so we don't have to do it for every iteration below
    one_minus_k = one(V) - k

    d = Dict(
        fe => sum(p[1].second * p[2].second for p in ps if (p[1].first ∩ p[2].first) == fe; init=zero(V)) / one_minus_k
        for fe in focal_elements
    )

    Z = BPA(d)
    setconflict!(Z, k)

    # Due to rounding, the total mass might be slightly greater than one after combination

    return Z
end

"""
    combine_yager(X::BPA, Y::BPA; conflict=true)

Combine two BPAs using Yager's Rule of Combination.

See also: [`combine_dempster`](@ref), [`bpa`](@ref).
"""
function combine_yager(X::BPA{K,V}, Y::BPA{K,V}; conflict=true) where {K,V}
    # Calculate the cross product of both mass assignments
    ps = Iterators.product(X, Y)

    # Get all focal elements
    focal_elements = (focalelements(X) ∪ focalelements(Y))

    # Compute K, the mass of conflict
    k = sum(p[1].second * p[2].second for p in ps if isempty(p[1].first ∩ p[2].first))

    Z = BPA(
        (
            fe => sum(p[1].second * p[2].second for p in ps if (p[1].first ∩ p[2].first) == (fe ∩ fe); init=0)
            for fe in focal_elements
        )...
    )

    all_keys = reduce(union, keys(Z))

    if conflict
        Z[all_keys] += k
    end

    setconflict!(Z, k)

    return Z
end