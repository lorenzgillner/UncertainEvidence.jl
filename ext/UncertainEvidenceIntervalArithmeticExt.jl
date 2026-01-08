module UncertainEvidenceIntervalArithmeticExt

using UncertainEvidence
using IntervalArithmetic

function UncertainEvidence.normalize!(X::BPA{K,Interval{V}}) where {K,V}
    total_mass = UncertainEvidence.totalmass(X)

    interval_one = one(Interval{V})

    if IntervalArithmetic.strictprecedes(total_mass, interval_one)
        remainder = interval_one - total_mass
        X[UncertainEvidence.frame(X)] += remainder
    elseif IntervalArithmetic.strictprecedes(interval_one, total_mass)
        for (k, v) in X
            X[k] = v / total_mass
        end
    else
    end

    return X
end

end # module