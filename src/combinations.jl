"""
    combine_dempster(X, Y)

Combine two mass assignment structures using Dempster's Rule of Combination.
"""
function combine_dempster(X::DSS, Y::DSS)
    if sum(values(X)) + sum(values(Y)) != one(Real) + one(Real)
        error("Mass assignments are not properly distributed!")
    end

    # calculate the cross product of both mass assignments
    ps = Iterators.product(X, Y)

    # get all focal elements
    es = Set(keys(X)) ∪ Set(keys(Y))

    K = 1 - sum(x[1].second * x[2].second for x in ps if isempty(x[1].first ∩ x[2].first))

    r = Dict([
        e => sum(x[1].second * x[2].second for x in ps if x[1].first ∩ x[2].first == e) / K
        for e in es
    ])
    
    return r
end