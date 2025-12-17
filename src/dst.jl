const BaseType = Number;

wrapset(x::T) where {T<:AbstractArray} = Set(x);
wrapset(x) = x isa AbstractSet ? x : Set([x]);

"""
    BPA{K,V} where {K<:Any, V<:Number}

A basic probability assessment (BPA) is the foundational data structure
for calculations in the context of the Dempster-Shafer theory (DST).

Use `bpa` to create automatically normalized BPAs.
Otherwise, use `redistribute!` for normalization.

See also: [`bpa`](@ref), [`redistribute!`](@ref).
"""

# BPA is a subtype of AbstractDict
struct BPA{K<:Any,V<:BaseType} <: AbstractDict{K,V}
    self::Dict{K,V}
    Ω::K

    function BPA{K,V}(d::Dict{K,V}) where {K<:AbstractSet,V<:BaseType}
        Ω = reduce(union, keys(d));
        new{K,V}(d, Ω);
    end

    function BPA{K,V}(d::Dict{K,V}) where {K<:Any,V<:BaseType}
        T = reduce(promote_type, [k isa AbstractArray ? eltype(k) : typeof(k) for k in keys(d)]);

        dd = Dict{Set{T},V}();

        for (k, v) in d
            dd[wrapset(k)] = v;
        end

        Ω = reduce(union, keys(dd));
        new{Set{T},V}(dd, Ω);
    end
end

# General construction
BPA() = BPA{Any,BaseType}(Dict{Any,BaseType}())
BPA(d::Dict{K,V}) where {K,V} = BPA{K,V}(d)
BPA(ps::Pair...) = BPA(ps)
BPA(ps::Pair{K,V}...) where {K,V} = BPA(Dict{K,V}(ps))
BPA{K,V}(ps::Pair{K,V}...) where {K,V} = BPA{K,V}(Dict{K,V}(ps))
BPA(itr) = BPA(Dict(itr))
BPA{K,V}(itr) where {K,V} = BPA(Dict{K,V}(itr))

# AbstractDict interface
Base.length(X::BPA) = length(X.self)

Base.iterate(X::BPA) = iterate(X.self)
Base.iterate(X::BPA, i) = iterate(X.self, i)

Base.keys(X::BPA) = keys(X.self)
Base.values(X::BPA) = values(X.self)
Base.pairs(X::BPA) = pairs(X.self)

Base.getindex(X::BPA{Set{K},V}, k::K) where {K,V} = getindex(X.self, wrapset(k))
Base.getindex(X::BPA{Set{K},V}, k::Set{K}) where {K,V} = getindex(X.self, k)
Base.setindex!(X::BPA{K,V}, v::V, k::K) where {K,V} = (X.self[k] = v)

# BPA-specific accessor aliases
focalelements(X::BPA) = keys(X.self)
masses(X::BPA) = values(X.self)

# Display
Base.display(X::BPA{Set{K},V}) where {K,V} = begin
    println("BPA{$K, $V} with $(length(X)) entries:");
    for (k, v) in X
        if k != X.Ω
            println("    {$(join(k, ','))} => $v");
        end
    end
    println("    {$(join(X.Ω, ','))} => $(X[X.Ω])");
end

"""
    bpa(X...)

Create a normalized basic probability assignment (BPA) structure
from pairs of mass assignments `X`.

# Examples
```juliadoctest
julia> A = bpa("a" => 0.1, "b" => 0.2)
BPA{Char, Float64} with 3 entries:
    {'a'}      => 0.1
    {'b'}      => 0.2
    {'a', 'b'} => 0.7
```

See also: [`BPA`](@ref).
"""
function bpa(X...)
    return redistribute!(BPA(X...))
end

"""
    redistribute!(X)

Normalize a BPA so that the sum of all mass assignments is equal to 1.

See also: [`BPA`](@ref), [`bpa`](@ref).
"""
function redistribute!(X::BPA{K,V}) where {K,V}
    real_one = one(BaseType)

    Ω = reduce(union, focalelements(X))

    if Ω ∉ focalelements(X)
        X[Ω] = zero(V)
    end

    total_mass = sum(masses(X))

    if total_mass < real_one
        # If the sum of all focal elements, including Ω, is less
        # than 1, the remainder must be added to Ω.
        remainder = real_one - total_mass
        X[Ω] += remainder
    elseif total_mass > real_one
        # Normalize masses if their sum is greater than one;
        # in the case of intervals, strict relations apply
        for (k, v) in X
            X[k] = v / total_mass
        end
    else
        # `total_mass` is equal to 1; do nothing.
    end

    return X
end

"""
    bel(e, X)

Calculate the belief value for a focal element `e` in a BPA `X`.

See also: [`BPA`](@ref), [`pls`](@ref).
"""
function bel(e, X::BPA)
    rv = zero(BaseType)

    for (k, v) in X
        if issubset(k, e)
            rv += v
        end
    end

    return rv
end

"""
    pls(e, X)

Calculate plausibility value for a focal element `e` in a BPA `X`.

See also: [`BPA`](@ref), [`bel`](@ref).
"""
function pls(e, X::BPA)
    rv = zero(BaseType)

    for (k, v) in X
        if !isdisjoint(k, e)
            rv += v
        end
    end

    return rv
end