using StructuralIdentifiability, RationalFunctionFields

ode = @ODEmodel(
    x1'(t) = -(k21 + k31 + k01) * x1(t) + k12 * x2(t) + k13 * x3(t)  + u(t),
    x2'(t) = k21 * x1(t) - k12 * x2(t),
    x3'(t) = k31 * x1(t) - k13 * x3(t),
    y(t) = x1(t)
)

original_generators = StructuralIdentifiability.initial_identifiable_functions(ode, prob_threshold=0.99, with_states=false)[1][1]
@info original_generators
original_generators = [f // original_generators[1] for f in original_generators[2:end]]
println("Original generators:")
for f in original_generators
    println("\t$f")
end

simplified_generators = simplified_generating_set(RationalFunctionField(original_generators))
println("Simplified generators:")
for f in simplified_generators
    println("\t$f")
end

println("With minimization applied:")
simplified_generators = simplified_generating_set(RationalFunctionField(original_generators), enforce_minimality = true)
for f in simplified_generators
    println("\t$f")
end
