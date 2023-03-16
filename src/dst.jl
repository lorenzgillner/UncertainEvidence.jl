function redistribute(em::Dict{T,S}) where {S} where {T}
    Ω = reduce(∪, keys(em))
    if Ω ∉ keys(em)
        em[Ω] = 0.0
    end
    vs = sum(values(em))
    if vs < 1.0
        em[Ω] = 1.0 - vs
    elseif vs > 1.0 || 1.0 ∈ vs
        em = Dict(zip(keys(em), values(em) ./ vs))
    end
    return em
end

function combine_dempster(X::Dict, Y::Dict)
    X = redistribute(X)
    Y = redistribute(Y)
    ps = Iterators.product(X, Y)
    es = Set(keys(X)) ∪ Set(keys(Y))
    k = 1 - sum(x[1].second * x[2].second for x in ps if isempty(x[1].first ∩ x[2].first))
    r = Dict([
        e => sum(x[1].second * x[2].second for x in ps if x[1].first ∩ x[2].first == e) / k
        for e in es
    ])
    return r
end

"""
    bel(e, X)

Calculate belief of focal element `e` in set `X`.
"""
function bel(e, X::Dict)
    z = zero(eltype(X).types[2])
    for x in X
        if issubset(x.first, e)
            z += x.second
        end
    end
    return z
end

"""
    pls(e, X)

Calculate plausibility of focal element `e` in set `X`.
"""
function pls(e, X::Dict)
    z = zero(eltype(X).types[2])
    for x in X
        if !isdisjoint(x.first, e)
            z += x.second
        end
    end
    return z
end