import Pkg; Pkg.activate(joinpath(@__DIR__, "..", "env"))

using Revise
using StructuralIdentifiability
using RationalFunctionFields
using ParamPunPam
using Nemo
using Groebner
using BenchmarkTools
using PrettyTables
using DataStructures

skipped = ["NFkB", "Covid2", "Akt", "TumorPillis", "TumorHu", "LeukaemiaLeon2021", "MAPK_5out_bis", "cLV2", "QWWC"]

if !isdefined(Main, :PROBLEM_CACHE)
PROBLEM_CACHE = Dict{Any, Any}()
end

power_sum(x, k) = sum(x .^ k)
power_sums(x, k) = [power_sum(x, i) for i in 1:k]

function max_terms_of_coeffs(gb)
    isnothing(gb) && return nothing
    N, D = 0, 0
    for f in gb
        N, D = argmax(sum, vcat((N,D), map(c -> (length(numerator(c)), length(denominator(c))), collect(coefficients(f)))))
    end
    N, D
end

function max_degree_of_coeffs(gb)
    isnothing(gb) && return nothing
    N, D = 0, 0
    for f in gb
        N, D = argmax(sum, vcat((N,D), map(c -> (total_degree(numerator(c)), total_degree(denominator(c))), collect(coefficients(f)))))
    end
    N, D
end

begin
R, x = polynomial_ring(QQ, :x => (1:5,))
RR, y = polynomial_ring(fraction_field(R), :y => (1:5,))
gb = [y[1] - prod(x) // (sum(x) + x[1]^5), x[1]^8*y[2] - sum(x)^2, sum(x)*y[3] - 1 // sum(x)^2]
@assert max_terms_of_coeffs(gb) in ((length(sum(x)^2), 1), (1, length(sum(x)^2)))
@assert max_degree_of_coeffs(gb) == (5, 5)
end

function problems_simplification_identifiable_funcs(;with_states=true)
    # get `benchmarks`
    include(joinpath(dirname(dirname(pathof(StructuralIdentifiability))), "benchmarking", "benchmarks.jl"))
    global benchmarks
    odes = []
    
    for name in collect(keys(benchmarks))
        if any(s -> occursin(s, string(name)), skipped)
            @error "Skipping $name"
            continue
        end
        if !isempty(ARGS) && !occursin(ARGS[1], string(name))
            continue
        end
        push!(odes, benchmarks[name])
    end
    
    problems = []
    for ode in odes
        name = string(ode[:name], " with_states=$(with_states)")
        @info "" name
	if haskey(PROBLEM_CACHE, name)
	    push!(problems, PROBLEM_CACHE[name])
	    continue
    	end
	ode = ode[:ode]
        id_funcs =
            StructuralIdentifiability.initial_identifiable_functions(ode, prob_threshold=0.99, with_states=with_states)[1]
        problem = Dict(:name => name, :funcs => id_funcs)
        PROBLEM_CACHE[name] = problem
	push!(problems, problem)
    end
    problems
end

function problems_simplification_power_sums()
    problems = []
    for n in 4:7
        for k in (n - 1):(n + 1)
            R, x = polynomial_ring(Nemo.QQ, :x => (1:n,))
	    funcs = [vcat(R(1), power_sums(x, k))]
            name = "power_sums_$(n)_$(k)"
            problem = Dict(:name => name, :funcs => funcs)
            push!(problems, problem)
        end
    end
    problems
end

function a_single_run_of_implementation(funcs)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    D = (typemax(Int), typemax(Int))
    t = @elapsed RationalFunctionFields.groebner_basis_coeffs(rff, ordering=DegRevLex(), up_to_degree=D)
    @assert length(rff.mqs.cached_groebner_bases) == 1
    gb = first(values(rff.mqs.cached_groebner_bases))
    n_spec, n_reduce = rff.mqs.stats.n_spec_mod_p, rff.mqs.stats.n_red_mod_p

    @info "Time:" t
    @info "Number of evaluations:" n_spec n_reduce
    @info "Interpolated max. degree:" max_degree_of_coeffs(gb)
    @info "Interpolated max. number of terms:" max_terms_of_coeffs(gb)
    return Dict(
        :n_spec=>n_spec,
        :n_red=>n_reduce,
        :our_deg=>max_degree_of_coeffs(gb),
        :terms=>max_terms_of_coeffs(gb),
        :time=>t
    )
end

function compute_max_total_degrees(funcs)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    degs = paramgb_only_degrees(rff.mqs, ordering=DegRevLex())
    N, D = 0, 0
    for i in degs
	N, D = argmax(sum, vcat((N,D), i))
    end
    return Dict{Any,Any}(
        :max_deg => (N, D)
    )
end

function try_to_compute_full_gb(funcs)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    t = @elapsed gb = paramgb(rff.mqs)
    n_spec, n_reduce = rff.mqs.stats.n_spec_mod_p, rff.mqs.stats.n_red_mod_p
    return Dict(
        :n_spec_full=>n_spec,
        :n_red_full=>n_reduce,
        :terms_full=>max_terms_of_coeffs(gb),
        :time_full=>t
    )
end

function compute_sharp_n_spec(funcs; deg=nothing)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    t = @elapsed gb = paramgb(rff.mqs, up_to_degree=deg)
    n_spec, n_reduce = rff.mqs.stats.n_spec_mod_p, rff.mqs.stats.n_red_mod_p
    return Dict(
        :n_spec_sharp=>n_spec
    )
end

function compute_time_per_eval(funcs)
    # p = Nemo.Native.GF(2^62 + 135)
    p = Nemo.Native.GF(2^30 + 3)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    StructuralIdentifiability.ParamPunPam.reduce_mod_p!(rff.mqs, p)
    point = rand(p, length(Nemo.gens(StructuralIdentifiability.ParamPunPam.parent_params(rff.mqs))))
    eqs = StructuralIdentifiability.ParamPunPam.specialize_mod_p(rff.mqs, point)
    print("groebner (rand.)   ")
    b1 = @belapsed Groebner.groebner($eqs, linalg=:randomized, ordering=Groebner.DegRevLex())
    display(b1)
    print("groebner (det.)    ")
    b2 = @belapsed Groebner.groebner($eqs, linalg=:deterministic, ordering=Groebner.DegRevLex())
    display(b2)
    print("groebner_apply!    ")
    trace, _ = Groebner.groebner_learn(eqs, ordering=DegRevLex())
    b3 = @belapsed Groebner.groebner_apply!($trace, $eqs, ordering=Groebner.DegRevLex())
    display(b3)
    return Dict{Any,Any}(
        :time_per_gb => b1,
    	:time_per_gb_det => b2,
        :time_per_apply => b3
    )
end

function basic_stats(funcs)
    n = length(gens(parent(funcs[1][1])))
    N, D = 0, 0
    for arr in funcs
	den = arr[1]
	degs = map(num -> (total_degree(num), total_degree(den)), arr[2:end])
	N, D = argmax(sum, vcat((N,D), degs))
    end
    return Dict{Any,Any}(
        :input_n => n,
        :input_deg => (N, D)
    )
end
 
function push_to_table!(table, res)
    for k in keys(table)
        if k in keys(res)
            push!(table[k], res[k])
        else
            push!(table[k], nothing)
        end
    end
end

problems = vcat(
    # problems_simplification_power_sums(),
    problems_simplification_identifiable_funcs(with_states=false),
)

table = OrderedDict(:name => [], :input_n => [], :input_deg => [], :n_spec => [], :n_red => [], :our_deg => [], :max_deg => [], :terms => [], :time => [], :time_per_gb => [], :time_per_gb_det => [], :time_per_apply => [], :terms_full => [], :n_spec_full => [], :n_red_full => [], :time_full => [], :n_spec_sharp => [])

for problem in problems
    name = problem[:name]
    funcs = problem[:funcs]

    println("============= $name ================")
    res0 = basic_stats(funcs)
    res1 = a_single_run_of_implementation(funcs)
    res2 = compute_max_total_degrees(funcs)
    res3 = compute_time_per_eval(funcs)
    # res4 = try_to_compute_full_gb(funcs)
    res4 = Dict()
    res5 = compute_sharp_n_spec(funcs, deg = max.(1, res1[:our_deg]))
    res = merge(res0, res1, res2, res3, res4, res5)
    res[:name] = name
    push_to_table!(table, res)
    println("====================================")
end

pretty_table(reduce(hcat, [v for (k, v) in table]), column_labels=string.(collect(keys(table))))

#=
begin
table[:name] = map(c -> chopsuffix(c, " with_states=false"), table[:name])
table_6 = OrderedDict(k => table[k] for k in [:name, :input_n, :input_deg, :max_deg, :our_deg, :terms, :n_red, :n_spec, :n_spec_sharp, :time]);
perm = sort(collect(1:length(table_6[:name])), by=i -> table_6[:name][i]);
for key in collect(keys(table_6))
   table_6[key] = table_6[key][perm]
end
table_6[:name] = map(c -> length(c) > 40 ? c[1:40] : c, table_6[:name])
pretty_table(reduce(hcat, [v for (k, v) in table_6]), column_labels=string.(collect(keys(table_6))), fit_table_in_display_vertically=false)
pretty_table(reduce(hcat, [v for (k, v) in table_6]), column_labels=string.(collect(keys(table_6))), backend=:markdown)
end
=#


=#
