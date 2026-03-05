using StructuralIdentifiability, RationalFunctionFields

ode = @ODEmodel(
    S'(t) = -beta * S(t) * I(t) / N - delta * beta * S(t) * T(t) / N,
    I'(t) =
        beta * S(t) * I(t) / N + delta * beta * S(t) * T(t) / N - (alpha + gamma) * I(t),
    T'(t) = gamma * I(t) - nu * T(t),
    y1(t) = T(t)
)

original_generators = StructuralIdentifiability.initial_identifiable_functions(
    ode,
    prob_threshold = 0.99,
    with_states = false,
)[1][1]
original_generators = [f // original_generators[1] for f in original_generators[2:end]]
println("Original generators:")
for f in original_generators
    println("\t$f")
end

simplified_generators =
    simplified_generating_set(RationalFunctionField(original_generators))
println("Simplified generators:")
for f in simplified_generators
    println("\t$f")
end
