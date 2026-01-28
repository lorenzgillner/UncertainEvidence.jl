"""
    bel(A, X::BPA)

Calculate the belief value for a focal element `A` in a BPA `X`.

See also: [`BPA`](@ref), [`pls`](@ref).
"""
function bel(A, X::BPA{K,V}) where {K,V}
    belief = zero(V)

    for (B, MB) in X
        if issubset(B, A)
            belief += MB
        end
    end

    return belief
end

"""
    pls(A, X::BPA)

Calculate plausibility value for a focal element `A` in a BPA `X`.

See also: [`BPA`](@ref), [`bel`](@ref).
"""
function pls(A, X::BPA{K,V}) where {K,V}
    plausibility = zero(V)

    for (B, MB) in X
        if !isdisjoint(B, A)
            plausibility += MB
        end
    end

    return plausibility
end

"""
    isnormal(X::BPA)

Check whether a BPA is normal, i.e. the total mass is equal to one.

See also: [`normalize!`](@ref), [`BPA`](@ref).
"""
isnormal(X::BPA{K,V}) where {K,V} = totalmass(X) ≈ one(V) # XXX Danger!

"""
    normalize!(X::BPA)

Normalize a BPA so that the sum of all mass assignments is equal to one.

See also: [`isnormal`](@ref), [`BPA`](@ref).
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