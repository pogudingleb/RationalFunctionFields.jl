using StructuralIdentifiability, RationalFunctionFields

ode = @ODEmodel(
    S'(t) = -beta * I(t) * S(t) / N - u(t) * S(t) / N,
    L'(t) = beta * I(t) * S(t) / N - alpha * L(t),
    I'(t) = alpha * L(t) - gamma * I(t) + sigma * Q(t),
    Q'(t) = (1 - eta) * gamma * I(t) - sigma * Q(t),
    y(t) = I(t) / N
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
