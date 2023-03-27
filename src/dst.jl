const DSS{T,S} = Dict{T,S} where {T<:Any, S<:Real}

DSS() = DSS{Any, Real}()
DSS(ps::Pair{T,S}...) where {T,S} = DSS{T,S}(ps)
function DSS(it)#::Iterators.Zip)
    temp = DSS()
    sizehint!(temp, sizeof(DSS)*length(it))
    for (k,v) in it
        temp[k] = v
    end
    return temp
end

"""
    dss(x...)

Create a normalized Dempster-Shafer structure (DSS) from pairs of focal elements and mass assignments.

# Examples
```juliadoctest
julia> A = dss(Set("a") => 0.1, Set("b") => 0.2)
Dict{Set{Char}, Float64} with 3 entries:
  Set(['a'])      => 0.1
  Set(['b'])      => 0.2
  Set(['a', 'b']) => 0.7
```
"""
function dss(X...)
    if Any in eltype(X).types
        @warn "Focal elements of type `Any` detected! Set operations may fail."
    end
    redistribute!(DSS(X...))
end

"""
    redistribute(X)

Normalize a Dempster-Shafer structure so that the sum of all mass assignments equals 1.
"""
function redistribute!(X::DSS)
    Ω = reduce(∪, keys(X))
    
    if Ω ∉ keys(X)
        X[Ω] = zero(Real)
    end

    vs = sum(values(X))

    if vs < one(Real)
        # If the sum of all focal elements, including Ω, is less
        # than 1, the remainder must be added to Ω.
        X[Ω] += one(Real) - vs
    elseif vs > one(Real) || one(Real) ∈ vs
        # Normalize masses if their sum is greater than one;
        # in some cases (e.g. intervals) the containment of 1
        # is the sufficient criterion for normalization.
        for (k,v) in X
            X[k] = v / vs
        end
    else
        # `vs` must be equal to 1; do nothing.
    end

    return X
end

"""
    bel(e, X)

Calculate the belief value for a focal element `e` in a Dempster-Shafer structure `X`.
"""
function bel(e, X::DSS)
    z = zero(Real)

    for x in X
        if issubset(x.first, e)
            z += x.second
        end
    end

    return z
end

"""
    pls(e, X)

Calculate plausibility value for a focal element `e` in a Dempster-Shafer structure `X`.
"""
function pls(e, X::DSS)
    z = zero(Real)

    for x in X
        if !isdisjoint(x.first, e)
            z += x.second
        end
    end

    return z
end