# Contains example for Section 4.2

using StructuralIdentifiability, RationalFunctionFields, ParamPunPam, Nemo

# 0. Get the defintion of the example  
include(
    joinpath(
        dirname(dirname(pathof(StructuralIdentifiability))),
        "benchmarking",
        "benchmarks.jl",
    ),
)
model = benchmarks[:SEIR34]

# 1. Get the original set of generators
funcs = StructuralIdentifiability.initial_identifiable_functions(
    model[:ode],
    prob_threshold = 0.99,
    with_states = false,
)[1]
println("=======\nSection 4.2\n=======\n")
println(
    "1. Original generators:\n",
    RationalFunctionFields.dennums_to_fractions(funcs),
    "\n",
)

# 2. Compute the coefficients of GB up to degree d=2 
#    (or, (d_num,d_den) = (1,1))
rff = RationalFunctionField(funcs)
deg = (1, 1)
gb = ParamPunPam.paramgb(rff.oms, up_to_degree = deg)
cfs = reduce(vcat, map(f -> collect(coefficients(f)), gb))
cfs_are_enough =
    fields_equal(RationalFunctionField(cfs), RationalFunctionField(funcs), 0.99)
println("2. Coefficients up to degree d=$deg:\n", cfs)
println("2. They generate the original subfield: ", cfs_are_enough)

# 3. Compute the coefficients of GB up to degree d=4
#    (or, (d_num,d_den) = (2,2))
rff = RationalFunctionField(funcs)
deg = (2, 2)
gb = ParamPunPam.paramgb(rff.oms, up_to_degree = deg)
cfs = reduce(vcat, map(f -> collect(coefficients(f)), gb))
cfs_are_enough =
    fields_equal(RationalFunctionField(cfs), RationalFunctionField(funcs), 0.99)
println("\n3. Coefficients up to degree d=$deg:\n", cfs)
println("3. They generate the original subfield: ", cfs_are_enough)

polys = RationalFunctionFields.polynomial_generators(rff, 2)
println("\n4. Polynomial generators:\n", polys)

combined = vcat(reduce(vcat, cfs), polys)
combined = filter(!RationalFunctionFields.is_rational_func_const, combined)
println("\n5. Combined:\n", combined)

sorted = sort(unique(combined), lt = RationalFunctionFields.rational_function_cmp)
println("\n6. Sorted:\n", sorted)

simplified = RationalFunctionFields.beautiful_generators(RationalFunctionField(combined))
println("\n7. Filtered:\n", simplified)

println("\nFinal simplified generators:\n", join(string.(simplified), ", "))
