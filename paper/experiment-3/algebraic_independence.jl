using Nemo, LinearAlgebra

function are_algebraically_independent(ratfuncs; p=0.99)
    base_vars = gens(base_ring(parent(ratfuncs[1])))
    maxdeg = maximum([
        max(total_degree(numerator(f)), total_degree(denominator(f))) for
        f in ratfuncs
    ])
    # degree of the polynomial whose nonvanishing will be needed for correct result
    D = Int(ceil(2 * maxdeg * (length(ratfuncs) + 1)^3 * length(ratfuncs) / (1 - p)))
    eval_point = [Nemo.QQ(rand(1:D)) for x in base_vars]

    # Filling the jacobain for generators
    S = matrix_space(Nemo.QQ, length(base_vars), length(ratfuncs))
    J = zero(S)
    for (i, f) in enumerate(ratfuncs)
        for (j, x) in enumerate(base_vars)
            J[j, i] = evaluate(derivative(f, x), eval_point)
        end
    end
    LinearAlgebra.rank(J) == length(ratfuncs)
end
