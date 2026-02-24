import Combinatorics: powerset
import PrettyTables: pretty_table

wrapset(x::T) where {T<:AbstractArray} = Set(x)
wrapset(x::Tuple) = Set(x)
wrapset(x) = x isa AbstractSet ? x : Set([x])

deducetype(d::Dict) = reduce(promote_type, [k isa AbstractArray ? eltype(k) : typeof(k) for k in keys(d)])

"""
    BPA{K,V} where {K<:AbstractSet,V<:Number}

A basic probability assignment (BPA) is the foundational data structure
for calculations in the context of the Dempster-Shafer theory (DST).
"""
mutable struct BPA{K<:Any,V<:Number} <: AbstractDict{Set{K},V}
    m::Dict{Set{K},V} # TODO Use two arrays instead
    Ω::Set{K} # TODO Don't store the key twice; use a reference
    conflict::V

    function BPA{K,V}(d::Dict{Set{K},V}, Ω::Set{K}) where {K<:Any,V<:Number}
        new{K,V}(d, Ω, zero(V))
    end
end

function BPA(d::Dict{Set{K},V}, Ω::Set{K}) where {K<:Any,V<:Number}
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
        d[Ω] = total_mass < one(V) ? one(V) - total_mass : zero(V) # TODO handle intervals
    end

    BPA{K,V}(d, Ω)
end

# Constructor for BPA with optional Ω; deduce Ω from keys if not supplied
function BPA(d::Dict{Set{K},V}; Ω=missing) where {K<:Any,V<:Number}
    ω = ismissing(Ω) ? reduce(union, keys(d)) : Ω
    BPA(d, ω)
end

# Constructor for BPA from non-set keys
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

# Constructor for BPA from non-set keys with optional Ω
function BPA(d::Dict{K,V}; Ω=missing) where {K<:Any,V<:Number}
    ω = ismissing(Ω) ? reduce(union, wrapset.(keys(d))) : Ω
    BPA(d, ω)
end

BPA(d::Dict{K,V}, Ω::AbstractArray{K}) where {K<:Any,V<:Number} = BPA(d, wrapset(Ω))

BPA() = throw(ArgumentError("Can't create a BPA from nothing; the mass of ∅ must zero"))

# BPA(d::Dict; Ω=missing) = ismissing(Ω) ? BPA(d) : BPA(d, Ω)
BPA(ps::Pair...; Ω=missing) = BPA(Dict(ps); Ω=Ω)
BPA(k::Array{Set{K}}, v::Array{V}) where {K,V} = BPA(Pairs(k, v))

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

function setconflict!(X::BPA{K,V}, c::V) where {K,V}
    X.conflict = c
end

# Display function
Base.print(io::IO, X::BPA{K,V}; showzero=false, showconflict=true) where {K,V} = begin
    println(io, "BPA{$K, $V} with $(length(X)) entries:")
    header = ["Focal element", "Mass"]
    tabular = vcat((["$(join(k, ','))" v] for (k, v) in X if (!iszero(v) || showzero))...)
    tabular = tabular[sortperm(tabular[:,1], by=length),:]
    pretty_table(io, tabular; column_labels=header, alignment=[:l, :r], compact_printing=true)
    if showconflict
        println(io, "Conflict = $(X.conflict)")
    end
end