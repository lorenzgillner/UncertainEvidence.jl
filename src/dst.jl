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