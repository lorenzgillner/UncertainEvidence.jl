function distribute(em::Dict{T,S}) where {S} where {T}
    Ω = reduce(∪, keys(em))
    # Ω = T(reduce(∪, keys(em)))
    if Ω ∉ keys(em)
        em[Ω] = 0.0
    end
    vs = sum(values(em))
    if vs < 1.0
        em[Ω] = 1.0 - vs
    elseif vs > 1.0
        em = Dict(zip(keys(em), values(em) ./ vs))
    end
    return em
end

# function dss(em::Dict)
#     return distribute(em)
# end

# TODO other types

function combine_dempster(X::Dict, Y::Dict)
    X = distribute(X)
    Y = distribute(Y)
    ps = Iterators.product(X, Y)
    es = Set(keys(X)) ∪ Set(keys(Y))
    k = 1 - sum(x[1].second * x[2].second for x in ps if isempty(x[1].first ∩ x[2].first))
    Dict([
        e => sum(x[1].second * x[2].second for x in ps if x[1].first ∩ x[2].first == e) / k
        for e in es
    ])
end

function bel(e, X::Dict)
    sum([x.second for x in X if x.first ⊆ e])
end

function pls(e, X::Dict)
    sum([x.second for x in X if !isempty(x.first ∩ e)])
end
