import Combinatorics: powerset

wrapset(x::T) where {T<:AbstractArray} = Set(x)
wrapset(x) = x isa AbstractSet ? x : Set([x])

deducetype(d::Dict) = reduce(promote_type, [k isa AbstractArray ? eltype(k) : typeof(k) for k in keys(d)])

"""
    BPA{K,V} where {K<:AbstractSet, V<:Number}

A basic probability assignment (BPA) is the foundational data structure
for calculations in the context of the Dempster-Shafer theory (DST).

See also: [`bpa`](@ref), [`redistribute!`](@ref).
"""
struct BPA{K<:Any,V<:Number} <: AbstractDict{Set{K},V}
    m::Dict{Set{K},V}
    Ω::Set{K}

    function BPA{K,V}(d::Dict{Set{K},V}, Ω::Set{K}) where {K<:Any,V<:Number}
        ks = keys(d)

        all_combinations = wrapset.(collect(powerset(collect(Ω), 1, length(Ω) - 1)))

        implicit_combinations = setdiff(all_combinations, ks)

        for k in implicit_combinations
            if isdisjoint(k, Ω)
                throw(ArgumentError("Focal element $k is not a subset of Ω"))
            end

            d[k] = zero(V)
        end

        if !in(Ω, ks)
            d[Ω] = one(V) - sum(values(d))
        end

        if sum(values(d)) > one(V)
            @warn "Sum of masses is greater than one. Consider redistributing"
        end

        new{K,V}(d, Ω)
    end
end

BPA(d::Dict{Set{K},V}, Ω::Set{K}) where {K<:Any,V<:Number} = BPA{K,V}(d, Ω)

# Constructor for BPA with optional Ω; deduce Ω from keys if not supplied
function BPA(d::Dict{Set{K},V}; Ω=nothing) where {K<:Any,V<:Number}
    ω = isnothing(Ω) ? reduce(union, keys(d)) : Ω
    BPA(d, ω)
end

# Constructor for non-set keys
function BPA(d::Dict{K,V}, Ω::Set{K}) where {K<:Any,V<:Number}
    S = Set{deducetype(d)}

    if S != typeof(Ω)
        throw(ArgumentError("Provided Ω type $(typeof(Ω)) does not match inferred key type $S"))
    end

    dd = Dict{S,V}()

    for (k, v) in d
        dd[wrapset(k)] = v
    end

    BPA(dd, Ω)
end

BPA(d::Dict{K,V}, Ω::AbstractArray{K}) where {K<:Any,V<:Number} = BPA(d, wrapset(Ω))

# Constructor for non-set keys without Ω; infer Ω from keys
function BPA(d::Dict{K,V}) where {K<:Any,V<:Number}
    S = Set{deducetype(d)}

    dd = Dict{S,V}()

    for (k, v) in d
        dd[wrapset(k)] = v
    end

    BPA(dd)
end

BPA() = throw(ArgumentError("Cannot create an empty BPA"))

# Convenience constructors
BPA(d::Dict; Ω=nothing) = isnothing(Ω) ? BPA(d) : BPA(d, Ω)
BPA(ps::Pair...; Ω=nothing) = isnothing(Ω) ? BPA(Dict(ps)) : BPA(Dict(ps), Ω)

# BPA(ps::Pair...) = BPA(Dict(ps))
# BPA(ps::Pair{K,V}...) where {K,V} = BPA(Dict{K,V}(ps))
# BPA{K,V}(ps::Pair{K,V}...) where {K,V} = BPA{K,V}(Dict{K,V}(ps))

# BPA(itr) = BPA(Dict(itr))
# BPA(itr, Ω) = BPA(Dict(itr), Ω)
# BPA{K,V}(itr) where {K,V} = BPA(Dict{K,V}(itr))

# AbstractDict interface
Base.length(X::BPA) = length(X.m)

Base.iterate(X::BPA) = iterate(X.m)
Base.iterate(X::BPA, i) = iterate(X.m, i)

Base.keys(X::BPA) = keys(X.m)
Base.values(X::BPA) = values(X.m)
Base.pairs(X::BPA) = pairs(X.m)

Base.getindex(X::BPA{K,V}, k::K) where {K,V} = getindex(X.m, wrapset(k))
Base.getindex(X::BPA{K,V}, k::Array{K}) where {K,V} = (k == []) ? zero(K) : getindex(X.m, wrapset(k))
Base.getindex(X::BPA{K,V}, k::Set{K}) where {K,V} = (k == Set{K}()) ? zero(K) : getindex(X.m, k)

function Base.setindex!(X::BPA{K,V}, v::V, k::U) where {K,V,U<:Union{K,Array{K},Set{K}}}
    if isdisjoint(wrapset(k), omega(X))
        throw(ArgumentError("$sk is not a focal element of Ω"))
    end
    X.m[k] = v
end

Base.eltype(X::BPA) = eltype(X.m)

# BPA-specific accessor aliases
# TODO rename these
focalelements(X::BPA) = keys(X.m)
masses(X::BPA) = values(X.m)
omega(X::BPA) = X.Ω

# Display function
Base.display(X::BPA{Set{K},V}) where {K,V} = begin
    println("BPA{$K, $V} with $(length(X)) entries:")
    for (k, v) in X
        if k != omega(X)
            println("  {$(join(k, ", "))} => $v")
        end
    end
    println("  {$(join(X.Ω, ", "))} => $(X[X.Ω])")
    # TODO Align lines properly
end

# TODO BPA for numeric sets as focal elements (see "earthquake example")

# """
#     bpa(X...)

# Create a normalized basic probability assignment (BPA) structure
# from pairs of mass assignments `X`.

# # Examples
# ```juliadoctest
# julia> A = bpa("a" => 0.1, "b" => 0.2)
# BPA{Char, Float64} with 3 entries:
#     {'a'}      => 0.1
#     {'b'}      => 0.2
#     {'a', 'b'} => 0.7
# ```

# See also: [`BPA`](@ref).
# """
# function bpa(X...)
#     return redistribute!(BPA(X...))
# end

"""
    redistribute!(X)

Normalize a BPA so that the sum of all mass assignments is equal to 1.

See also: [`BPA`](@ref), [`bpa`](@ref).
"""
function redistribute!(X::BPA{K,V}) where {K,V}
    real_one = one(Number)

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