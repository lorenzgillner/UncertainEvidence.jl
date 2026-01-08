"""
    bel(A, X)

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
    pls(A, X)

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