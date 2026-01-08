import Combinatorics: powerset

wrapset(x::T) where {T<:AbstractArray} = Set(x)
wrapset(x::Tuple) = Set(x)
wrapset(x) = x isa AbstractSet ? x : Set([x])

deducetype(d::Dict) = reduce(promote_type, [k isa AbstractArray ? eltype(k) : typeof(k) for k in keys(d)])

"""
    BPA{K,V} where {K<:AbstractSet,V<:Number}

A basic probability assignment (BPA) is the foundational data structure
for calculations in the context of the Dempster-Shafer theory (DST).

See also: [`bpa`](@ref), [`normalize!`](@ref).
"""
struct BPA{K<:Any,V<:Number} <: AbstractDict{Set{K},V}
    m::Dict{Set{K},V}
    Ω::Set{K} # TODO Don't store the key twice; use a reference to it in `m` instead

    function BPA{K,V}(d::Dict{Set{K},V}, Ω::Set{K}) where {K<:Any,V<:Number}
        all_combinations = wrapset.(collect(powerset(collect(Ω), 1, length(Ω) - 1)))

        implicit_combinations = setdiff(all_combinations, keys(d))

        for k in implicit_combinations
            if isdisjoint(k, Ω)
                throw(ArgumentError("Focal element $k is not a subset of Ω"))
            end

            d[k] = zero(V)
        end

        total_mass = sum(values(d))

        if !haskey(d, Ω)
            d[Ω] = total_mass < one(V) ? one(V) - total_mass : zero(V)
        end

        # if sum(values(d)) > one(V)
        #     @warn "Sum of masses is greater than one"
        # end

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

BPA() = throw(ArgumentError("Can't create a BPA from nothing; the mass of ∅ must zero"))

# Convenience constructors
BPA(d::Dict; Ω=nothing) = isnothing(Ω) ? BPA(d) : BPA(d, Ω)
BPA(ps::Pair...; Ω=nothing) = isnothing(Ω) ? BPA(Dict(ps)) : BPA(Dict(ps), Ω)

# TODO BPA type for sets of real numbers (see "earthquake example")

# AbstractDict interface
Base.length(X::BPA) = length(X.m)

Base.iterate(X::BPA) = iterate(X.m)
Base.iterate(X::BPA, i) = iterate(X.m, i)

Base.keys(X::BPA) = keys(X.m)
Base.values(X::BPA) = values(X.m)
Base.pairs(X::BPA) = pairs(X.m)

Base.getindex(X::BPA{K,V}) where {K,V} = zero(V)
Base.getindex(X::BPA{K,V}, k::K) where {K,V} = getindex(X.m, wrapset(k))
Base.getindex(X::BPA{K,V}, ks::K...) where {K,V} = getindex(X.m, wrapset(ks))
Base.getindex(X::BPA{K,V}, k::Array{K}) where {K,V} = (k == []) ? zero(V) : getindex(X.m, wrapset(k))
Base.getindex(X::BPA{K,V}, k::Set{K}) where {K,V} = (k == Set{K}()) ? zero(V) : getindex(X.m, k)

function Base.setindex!(X::BPA{K,V}, v::V, k::U) where {K,V,U<:Union{K,Array{K},Set{K},Tuple{K}}}
    if isdisjoint(wrapset(k), frame(X))
        throw(ArgumentError("$sk is not a valid focal element"))
    end
    X.m[wrapset(k)] = v
end

Base.eltype(X::BPA) = eltype(X.m)

# BPA-specific accessor aliases
focalelements(X::BPA) = keys(X.m)
masses(X::BPA) = values(X.m)
totalmass(X::BPA) = sum(values(X.m))
frame(X::BPA) = X.Ω

isnormal(X::BPA{K,V}) where {K,V} = totalmass(X) == one(V)

# Display function
# TODO Use PrettyTables.jl
Base.display(X::BPA{K,V}) where {K,V} = begin
    println("BPA{$K, $V} with $(length(X)) entries:")
    for (k, v) in X
        if k != frame(X)
            println("  {$(join(k, ", "))} => $v")
        end
    end
    println("  {$(join(X.Ω, ", "))} => $(X[X.Ω])")
    # TODO Align lines properly
end

"""
    normalize!(X)

Normalize a BPA so that the sum of all mass assignments is equal to 1.

See also: [`BPA`](@ref), [`bpa`](@ref).
"""
function normalize!(X::BPA{K,V}) where {K,V}
    total_mass = totalmass(X)

    if total_mass < one(V)
        # If the sum of all focal elements, including Ω, is less
        # than one, the remainder must be added to Ω.
        remainder = one(V) - total_mass
        X[frame(X)] += remainder
    elseif total_mass > one(V)
        # Normalize masses if their sum is strictly greater than one
        for (k, v) in X
            X[k] = v / total_mass
        end
    else
        # `total_mass` must be equal to 1; do nothing.
    end

    return X
end