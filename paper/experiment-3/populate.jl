using Revise
using StructuralIdentifiability, RationalFunctionFields, Nemo

# get `benchmarks`
include(joinpath(dirname(dirname(pathof(StructuralIdentifiability))), "benchmarking", "benchmarks.jl"))

systems = []

for name in keys(benchmarks)
    push!(systems, Dict(:name => name))
end

prefix = "results"

templates = Dict(
    "simplify" => 
Dict(
"" =>
"""
using StructuralIdentifiability, RationalFunctionFields, Nemo

# get `benchmarks`
include(joinpath(dirname(dirname(pathof(StructuralIdentifiability))), \"benchmarking\", \"benchmarks.jl\"))
    
name = :{{name}}
with_states = {{with_states}}
ode = benchmarks[name][:ode]

funcs = StructuralIdentifiability.initial_identifiable_functions(ode, prob_threshold=0.99, with_states=with_states)[1]
# println("\\ninitial funcs (format: [[den, num...]])\\n", funcs)

rff = RationalFunctionFields.RationalFunctionField(funcs);

t = @elapsed simple_funcs = simplified_generating_set(rff)
println("\\ntime (s)\\n>>>", t)
        
println("\\nsimple funcs\\n", simple_funcs);

degs = map(f -> (total_degree(numerator(f)), total_degree(denominator(f))), simple_funcs)
trms = map(f -> (length(numerator(f)), length(denominator(f))), simple_funcs)
println("max deg = ", argmax(sum, degs))
println("max terms = ", argmax(sum, trms))
println("min deg = ", argmin(sum, degs))
println("min terms = ", argmin(sum, trms))
"""
),

    "maple_simplify" => 
Dict(
"_generate_gens" =>
"""
using StructuralIdentifiability, RationalFunctionFields, Nemo

# get `benchmarks`
include(joinpath(dirname(dirname(pathof(StructuralIdentifiability))), \"benchmarking\", \"benchmarks.jl\"))
    
name = :{{name}}
with_states = {{with_states}}
ode = benchmarks[name][:ode]

funcs = StructuralIdentifiability.initial_identifiable_functions(ode, prob_threshold=0.99, with_states=with_states)[1]

filename = string(name, \"_gens.mpl\")
open(joinpath(@__DIR__, filename), \"w\") do io
    fracs = RationalFunctionFields.dennums_to_fractions(funcs)
    fracs_str = join(map(s -> replace(s, \"//\" => \"/\", \"gamma\" => \"gama\", \"I\" => \"_I\"), map(string, fracs)), \", \")
    println(io, \"gens := {\", fracs_str, \"}:\")
end
""",

"_run_simplify" =>
"""
read "{{dir}}/{{name}}_gens.mpl";
read "{{dir}}/../../AllIdentifiableFunctions/ComputeIdentifiableFunctions.mpl";
        
simple := FilterGenerators(FieldToIdeal(gens));

interface(screenwidth=infinity):
lprint(\">>>\", simple);

quit;
"""
),

    "input_stats" =>
Dict(
"_collect_input_stats" =>
"""
using StructuralIdentifiability, RationalFunctionFields, Nemo

# get `benchmarks`
include(joinpath(dirname(dirname(pathof(StructuralIdentifiability))), "benchmarking", "benchmarks.jl"))
    
name = :{{name}}
with_states = {{with_states}}
ode = benchmarks[name][:ode]

funcs = StructuralIdentifiability.initial_identifiable_functions(ode, prob_threshold=0.99, with_states=with_states)[1]

filename = string(name, "_gens_stats.txt")
open(joinpath(@__DIR__, filename), "w") do io
    println(io, "format [[den, nums]]")
    println(io, "length(string(...)) = ", length(string(funcs)))
    println(io, "length(...), map(length, ...) = ", length(funcs), ", ", map(length, funcs))
    
    fracs = RationalFunctionFields.dennums_to_fractions(funcs)
    println(io, "\n>>>\n\nformat [num/den] (filter out constants)")
    fracs = filter(f -> (total_degree(numerator(f)) + total_degree(denominator(f))) > 0, fracs)
    println(io, "length(string(...)) = ", length(string(fracs)))
    println(io, "length(...) = ", length(fracs))
    degs = map(f -> (total_degree(numerator(f)), total_degree(denominator(f))), fracs)
    trms = map(f -> (length(numerator(f)), length(denominator(f))), fracs)
    println(io, "max deg(...) = ", argmax(sum, degs))
    println(io, "max terms(...) = ", argmax(sum, trms))
    println(io, "min deg(...) = ", argmin(sum, degs))
    println(io, "min terms(...) = ", argmin(sum, trms))

    xs = gens(base_ring(parent(fracs[1])))
    min_gen_deg_per_var = map(sum, map(x -> argmin(sum, vcat((Inf, Inf), map(f -> (total_degree(numerator(f)), total_degree(denominator(f))), filter(f -> degree(numerator(f), x) > 0 || degree(denominator(f), x) > 0, fracs)))), xs))
    min_deg_per_var = map(sum, map(x -> argmin(sum, vcat((Inf, Inf), map(f -> (degree(numerator(f), x), degree(denominator(f), x)), filter(f -> degree(numerator(f), x) > 0 || degree(denominator(f), x) > 0, fracs)))), xs))
    println(io, "\n")
    println(io, "vars = ", "["*join(map(string, xs), ", ")*"]")
    println(io, "min gen deg per var = ", min_gen_deg_per_var)
    println(io, "min deg per var = ", min_deg_per_var)
    
    println(io, "\n")
    println(io, "all degs = ", degs)
    println(io, "all terms = ", trms)
end
""",
),

    "independence" => 
Dict(
"" =>
"""
using StructuralIdentifiability, RationalFunctionFields, Nemo

# get `benchmarks`
include(joinpath(dirname(dirname(pathof(StructuralIdentifiability))), \"benchmarking\", \"benchmarks.jl\"))

include(joinpath(@__DIR__, "..", "..", "algebraic_independence.jl"))

name = :{{name}}
with_states = {{with_states}}
ode = benchmarks[name][:ode]

funcs = StructuralIdentifiability.initial_identifiable_functions(ode, prob_threshold=0.99, with_states=with_states)[1]
# println("\\ninitial funcs (format: [[den, num...]])\\n", funcs)

rff = RationalFunctionFields.RationalFunctionField(funcs);

t = @elapsed simple_funcs = simplified_generating_set(rff)
println("\\ntime (s)\\n>>>", t)
        
println("\\nsimple funcs\\n", simple_funcs);

println("independent: ", are_algebraically_independent(simple_funcs))
"""
),
)

method_ext = Dict(
    "simplify" => Dict("" => ".jl"),
    "maple_simplify" => Dict("_generate_gens" => ".jl", "_run_simplify" => ".mpl"),
    "input_stats" => Dict("_collect_input_stats" => ".jl"),
    "independence" => Dict("" => ".jl")
)

function foo()
    !isdir(joinpath(@__DIR__, prefix)) && mkdir(joinpath(@__DIR__, prefix))
    for system in systems
        name = system[:name]
        for method in collect(templates)
            method_name, method_scripts = method
            !isdir(joinpath(@__DIR__, prefix, method_name)) && mkdir(joinpath(@__DIR__, prefix, method_name))
            subs = (
                "{{name}}"  => string(name),
                "{{with_states}}"  => "false",
                "{{dir}}" => joinpath(@__DIR__, prefix, method_name),
            )
            for (script_name, script_template) in method_scripts
                content = replace(script_template, subs...)
                @assert !occursin("{{", content)
                file_name = string(name, script_name, method_ext[method_name][script_name])
                file_path = joinpath(@__DIR__, prefix, method_name, file_name)
                open(file_path, "w") do io
                    println(io, content)
                end
            end
        end
    end
end

foo()

