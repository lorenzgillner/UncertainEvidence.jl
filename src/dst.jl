"""
    BPA = Dict{T,S} where {T<:Any, S<:Real}

A basic probability assessment (BPA) is the foundational data structure
for calculations in the context of the Dempster-Shafer theory (DST).

Use `bpa` to create automatically normalized BPAs.
Otherwise, use `redistribute!` for normalization.

See also: [`bpa`](@ref), [`redistribute!`](@ref).
"""
const BPA{T,S} = Dict{T,S} where {T<:Any,S<:Real}

# general case
BPA() = BPA{Any,Real}()

# BPAs are basically `Dict`s, because pairs are more intuitive
BPA(ps::Pair{T,S}...) where {T,S} = BPA{T,S}(ps)

# create BPAs from an iterable object
function BPA(it)
    temp = BPA()
    sizehint!(temp, sizeof(BPA) * length(it))
    for (k, v) in it
        temp[k] = v
    end
    return temp
end

"""
    bpa(X...)

Create a normalized basic probability assignment (BPA) structure
from pairs of mass assignments `X`.

# Examples
```juliadoctest
julia> A = bpa(Set("a") => 0.1, Set("b") => 0.2)
Dict{Set{Char}, Float64} with 3 entries:
    Set(['a'])      => 0.1
    Set(['b'])      => 0.2
    Set(['a', 'b']) => 0.7
```

See also: [`BPA`](@ref).
"""
function bpa(X...)
    if Any in eltype(X).types
        @warn "Focal elements of type `Any` detected! Set operations might fail."
    end
    return redistribute!(BPA(X...))
end

"""
    redistribute!(X)

Normalize a BPA so that the sum of all mass assignments is equal to 1.

See also: [`BPA`](@ref), [`bpa`](@ref).
"""
function redistribute!(X::BPA)
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
        for (k, v) in X
            X[k] = v / vs
        end
    else
        # `vs` must be equal to 1; do nothing.
    end

    return X
end

"""
    bel(e, X)

Calculate the belief value for a focal element `e` in a BPA `X`.

See also: [`BPA`](@ref).
"""
function bel(e, X::BPA)
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

Calculate plausibility value for a focal element `e` in a BPA `X`.

See also: [`BPA`](@ref).
"""
function pls(e, X::BPA)
    z = zero(Real)

    for x in X
        if !isdisjoint(x.first, e)
            z += x.second
        end
    end

    return z
end