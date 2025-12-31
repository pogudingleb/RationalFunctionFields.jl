using Pkg;
Pkg.activate(@__DIR__)

using Revise
using StructuralIdentifiability
using RationalFunctionFields
using ParamPunPam
using Nemo
using Groebner
using Chairmarks
using PrettyTables
using DataStructures

include("gizmos.jl")

if !isdefined(Main, :PROBLEM_CACHE)
PROBLEM_CACHE = Dict{Any, Any}()
end

skipped = ["NFkB", "Covid2", "Akt", "TumorPillis", "TumorHu", "LeukaemiaLeon2021", "MAPK_5out_bis", "cLV2", "QWWC"]

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
        for j in i
            N = max(N, j[1])
            D = max(D, j[2])
        end
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
    p = Nemo.Native.GF(2^62 + 135)
    rff = RationalFunctionFields.RationalFunctionField(funcs);
    StructuralIdentifiability.ParamPunPam.reduce_mod_p!(rff.mqs, p)
    point = rand(p, length(Nemo.gens(StructuralIdentifiability.ParamPunPam.parent_params(rff.mqs))))
    eqs = StructuralIdentifiability.ParamPunPam.specialize_mod_p(rff.mqs, point)
    print("groebner           ")
    b1 = @b Groebner.groebner($eqs, ordering=Groebner.DegRevLex())
    display(b1)
    print("groebner_apply!    ")
    trace, _ = Groebner.groebner_learn(eqs, ordering=DegRevLex())
    b2 = @b Groebner.groebner_apply!($trace, $eqs, ordering=Groebner.DegRevLex())
    display(b2)
    return Dict{Any,Any}(
        :time_per_gb => b1.time,
        :time_per_apply => b2.time
    )
end

function basic_stats(funcs)
    n = length(gens(parent(funcs[1][1])))
    N, D = 0, 0
    for arr in funcs
        D = max(D, total_degree(arr[1]))
        for f in arr[2:end]
            N = max(N, total_degree(f))
        end
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

table = OrderedDict(:name => [], :input_n => [], :input_deg => [], :n_spec => [], :n_red => [], :our_deg => [], :max_deg => [], :terms => [], :time => [], :time_per_gb => [], :time_per_apply => [], :terms_full => [], :n_spec_full => [], :n_red_full => [], :time_full => [], :n_spec_sharp => [])

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
DEFAULT RUN

┌─────────────────────────────────────────────┬─────────┬───────────┬────────┬───────┬─────────┬──────────┬──────────┬───────────┬─────────────┬────────────────┐
│                                        name │ input_n │ input_deg │ n_spec │ n_red │ our_deg │  max_deg │    terms │      time │ time_per_gb │ time_per_apply │
├─────────────────────────────────────────────┼─────────┼───────────┼────────┼───────┼─────────┼──────────┼──────────┼───────────┼─────────────┼────────────────┤
│              SIWR original with_states=true │      11 │   (10, 0) │     18 │     1 │  (1, 0) │   (1, 0) │   (1, 1) │   2.40801 │   0.0008305 │        0.00043 │
│    Modified LV for testing with_states=true │       6 │    (5, 0) │     30 │     1 │  (2, 0) │   (2, 0) │   (2, 1) │  0.005589 │     6.23e-5 │        3.24e-5 │
│         Goodwin oscillator with_states=true │      11 │   (12, 7) │     62 │     1 │  (2, 2) │   (3, 3) │   (3, 2) │  0.453688 │   0.0004917 │      0.0001526 │
│                        HIV with_states=true │      15 │   (10, 1) │     58 │     1 │  (3, 3) │   (3, 3) │   (1, 1) │  0.692269 │   0.0004114 │        9.87e-5 │
│                SIRS forced with_states=true │      11 │   (17, 0) │     26 │     1 │  (2, 2) │   (2, 2) │   (1, 1) │ 0.0334116 │   0.0011353 │      0.0003902 │
│                Akt pathway with_states=true │      25 │    (7, 1) │     62 │     1 │  (2, 2) │   (2, 2) │   (3, 1) │   5.79549 │   0.0004707 │       0.000152 │
│  Chemical reaction network with_states=true │      12 │    (9, 0) │     18 │     1 │  (1, 0) │   (1, 0) │   (1, 1) │ 0.0222071 │   0.0004292 │      0.0001125 │
│                      SLIQR with_states=true │      10 │   (12, 0) │   1228 │     1 │  (7, 6) │   (7, 6) │ (23, 11) │  0.783162 │   0.0006723 │      0.0001341 │
│                         St with_states=true │      12 │    (7, 2) │    772 │     1 │  (4, 4) │ (11, 11) │ (29, 13) │   11.4842 │   0.0477545 │      0.0094507 │
│              Bilirubin2_io with_states=true │      11 │    (5, 0) │    644 │     1 │  (4, 2) │   (4, 2) │  (18, 6) │  0.904304 │   0.0018434 │      0.0004575 │
│             SIWR original with_states=false │       7 │  (24, 11) │     24 │     1 │  (1, 2) │   (1, 2) │   (1, 1) │   7.49199 │    0.442858 │      0.0696024 │
│   Modified LV for testing with_states=false │       4 │    (2, 3) │     30 │     1 │  (2, 0) │   (2, 0) │   (2, 1) │ 0.0062482 │     4.66e-5 │        2.71e-5 │
│        Goodwin oscillator with_states=false │       7 │   (14, 6) │     34 │     1 │  (2, 1) │   (2, 1) │   (2, 1) │ 0.0325746 │    0.000645 │      0.0001935 │
│                       HIV with_states=false │      10 │    (8, 3) │     56 │     1 │  (3, 2) │   (3, 2) │   (1, 1) │ 0.0604035 │     0.00062 │      0.0001358 │
│               SIRS forced with_states=false │       6 │   (15, 5) │     22 │     1 │  (2, 0) │   (2, 0) │   (1, 1) │   1.48899 │   0.0650385 │      0.0125389 │
│               Akt pathway with_states=false │      16 │  (15, 12) │     46 │     1 │  (1, 1) │   (4, 7) │   (3, 1) │  0.141464 │   0.0018867 │      0.0005381 │
│ Chemical reaction network with_states=false │       6 │    (5, 4) │     34 │     1 │  (1, 2) │   (1, 2) │   (1, 2) │ 0.0314471 │   0.0003445 │      0.0001182 │
│                     SLIQR with_states=false │       6 │    (6, 4) │    716 │     1 │  (7, 6) │   (7, 6) │ (14, 11) │  0.207696 │   0.0003166 │      0.0001008 │
│                        St with_states=false │       9 │  (18, 12) │    772 │     1 │  (4, 4) │  (11, 9) │  (26, 8) │   33.2523 │    0.336453 │      0.0123608 │
│             Bilirubin2_io with_states=false │       7 │    (4, 3) │    316 │     1 │  (4, 2) │   (4, 2) │  (12, 6) │  0.186747 │   0.0003214 │      0.0001259 │
└─────────────────────────────────────────────┴─────────┴───────────┴────────┴───────┴─────────┴──────────┴──────────┴───────────┴─────────────┴────────────────┘


|                            **name** | **our\_deg** | **max\_deg** |
|------------------------------------:|-------------:|-------------:|
|               SIWR original wt=true |       (1, 0) |       (1, 0) |
|     Modified LV for testing wt=true |       (2, 0) |       (2, 0) |
|                       Pharm wt=true |       (1, 2) |       (1, 2) |
|         SEAIJRC Covid model wt=true |       (2, 2) |       (2, 2) |
|      MAPK model (5 outputs) wt=true |       (1, 0) |       (1, 0) |
|      MAPK model (6 outputs) wt=true |       (1, 0) |       (1, 0) |
|          Goodwin oscillator wt=true |       (2, 2) |       (3, 3) |
|                         HIV wt=true |       (3, 3) |       (3, 3) |
|                 SIRS forced wt=true |       (2, 2) |       (2, 2) |
|                 Akt pathway wt=true |       (2, 2) |       (2, 2) |
|   Chemical reaction network wt=true |       (1, 0) |       (1, 0) |
|  CD8 T cell differentiation wt=true |       (2, 2) |       (2, 2) |
|                 LLW1987\_io wt=true |       (2, 2) |       (3, 2) |
|                    HIV2\_io wt=true |       (4, 4) |     (15, 14) |
|               Treatment\_io wt=true |       (2, 2) |       (3, 4) |
|        Biohydrogenation\_io wt=true |       (2, 1) |       (2, 4) |
|                       SLIQR wt=true |       (7, 6) |       (7, 6) |
|                          St wt=true |       (4, 4) |     (11, 11) |
|              Bilirubin2\_io wt=true |       (4, 2) |       (4, 2) |
|              SIWR original wt=false |       (1, 2) |       (1, 2) |
|    Modified LV for testing wt=false |       (2, 0) |       (2, 0) |
|                      Pharm wt=false |       (1, 2) |       (1, 2) |
|        SEAIJRC Covid model wt=false |       (2, 1) |       (2, 3) |
|     MAPK model (5 outputs) wt=false |       (1, 0) |      (1, 12) |
|     MAPK model (6 outputs) wt=false |       (1, 0) |       (1, 9) |
|         Goodwin oscillator wt=false |       (2, 1) |       (2, 1) |
|                        HIV wt=false |       (3, 2) |       (3, 2) |
|                SIRS forced wt=false |       (2, 0) |       (2, 0) |
|                Akt pathway wt=false |       (1, 1) |       (4, 7) |
|  Chemical reaction network wt=false |       (1, 2) |       (1, 2) |
| CD8 T cell differentiation wt=false |       (1, 1) |       (2, 3) |
|                LLW1987\_io wt=false |       (2, 0) |       (2, 0) |
|                   HIV2\_io wt=false |       (4, 2) |       (4, 2) |
|              Treatment\_io wt=false |       (2, 1) |       (4, 4) |
|       Biohydrogenation\_io wt=false |       (2, 2) |       (2, 2) |
|                      SLIQR wt=false |       (7, 6) |       (7, 6) |
|                         St wt=false |       (4, 4) |      (11, 9) |
|             Bilirubin2\_io wt=false |       (4, 2) |       (4, 2) |

julia> pretty_table(reduce(hcat, [v for (k, v) in table]), column_labels=string.(collect(keys(table))), fit_table_in_display_vertically=false)
┌──────────────────────┬─────────┬───────────┬────────┬───────┬─────────┬────────────┬──────────┬────────────┬─────────────┬────────────────┬────────────┬─────────────┬────────────┬───────────┬──────────────┐
│                 name │ input_n │ input_deg │ n_spec │ n_red │ our_deg │    max_deg │    terms │       time │ time_per_gb │ time_per_apply │ terms_full │ n_spec_full │ n_red_full │ time_full │ n_spec_sharp │
├──────────────────────┼─────────┼───────────┼────────┼───────┼─────────┼────────────┼──────────┼────────────┼─────────────┼────────────────┼────────────┼─────────────┼────────────┼───────────┼──────────────┤
│ p53 with_states=fals │      23 │    (9, 6) │     82 │     1 │  (4, 4) │     (4, 4) │   (1, 2) │    12.3169 │  0.00014903 │      4.8964e-5 │    nothing │     nothing │    nothing │   nothing │           62 │
│ HighDimNonLin with_s │      22 │    (1, 2) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │    1.31029 │   8.9872e-5 │      4.1492e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ SEIR 34 with_states= │       7 │    (5, 4) │     34 │     1 │  (2, 1) │     (2, 3) │   (2, 1) │    2.35504 │   5.1715e-5 │      2.5496e-5 │    nothing │     nothing │    nothing │   nothing │           34 │
│ Immune response to i │       8 │    (4, 1) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │ 0.00636277 │ 0.000113755 │      3.0333e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ Immune response to i │       7 │    (3, 1) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │    1.35283 │   5.5374e-5 │      2.0219e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ SIR 19 with_states=f │       6 │    (3, 3) │     20 │     1 │  (1, 1) │     (1, 1) │   (1, 1) │  0.0026336 │    3.353e-5 │      2.0461e-5 │    nothing │     nothing │    nothing │   nothing │           18 │
│ Goodwin oscillator w │       7 │   (14, 6) │     34 │     1 │  (2, 1) │     (2, 1) │   (2, 1) │  0.0128422 │ 0.000316511 │      8.7816e-5 │    nothing │     nothing │    nothing │   nothing │           34 │
│ CD8 T cell different │      13 │    (4, 4) │     22 │     1 │  (1, 1) │     (2, 3) │   (1, 1) │   0.174085 │ 0.000112436 │      4.5095e-5 │    nothing │     nothing │    nothing │   nothing │           18 │
│ Bruno2016 with_state │       6 │    (2, 1) │     24 │     1 │  (1, 0) │     (1, 0) │   (2, 1) │   0.144872 │    2.218e-5 │      1.2347e-5 │    nothing │     nothing │    nothing │   nothing │           22 │
│ Covid model (Chitnis │       8 │    (4, 5) │     38 │     1 │  (2, 2) │     (2, 3) │   (2, 1) │   0.174208 │   6.2863e-5 │      2.6941e-5 │    nothing │     nothing │    nothing │   nothing │           38 │
│ SIR 21 with_states=f │       6 │    (3, 2) │     20 │     1 │  (1, 1) │     (1, 1) │   (1, 1) │ 0.00260027 │   3.3709e-5 │      1.9403e-5 │    nothing │     nothing │    nothing │   nothing │           18 │
│ SEIRT with_states=fa │       4 │    (3, 2) │     34 │     1 │  (2, 1) │     (2, 1) │   (2, 1) │   0.143695 │   2.8734e-5 │      1.5868e-5 │    nothing │     nothing │    nothing │   nothing │           34 │
│ SEAIJRC Covid model  │       7 │  (31, 15) │     34 │     1 │  (2, 1) │     (2, 3) │   (2, 1) │    17.7058 │    0.763606 │       0.334579 │    nothing │     nothing │    nothing │   nothing │           34 │
│ MAPK model (6 output │      22 │    (9, 8) │     20 │     1 │  (1, 0) │     (1, 9) │   (1, 1) │   0.382679 │  0.00292888 │    0.000513414 │    nothing │     nothing │    nothing │   nothing │           16 │
│ Bilirubin2_io with_s │       7 │    (4, 3) │    316 │     1 │  (4, 2) │     (4, 2) │  (12, 6) │   0.067053 │ 0.000103133 │      4.1741e-5 │    nothing │     nothing │    nothing │   nothing │          278 │
│ SEIR_1_io with_state │       4 │    (7, 5) │     62 │     1 │  (2, 2) │     (3, 2) │   (4, 2) │  0.0062147 │    7.815e-5 │      3.1531e-5 │    nothing │     nothing │    nothing │   nothing │           62 │
│ EIHRD epidemiologica │      10 │  (22, 14) │    420 │     1 │  (4, 4) │    (10, 9) │  (14, 8) │    68.4431 │    0.671201 │       0.182964 │    nothing │     nothing │    nothing │   nothing │          342 │
│ SIRS forced with_sta │       6 │   (15, 5) │     22 │     1 │  (2, 0) │     (2, 0) │   (1, 1) │   0.449433 │   0.0281319 │     0.00686877 │    nothing │     nothing │    nothing │   nothing │           22 │
│ Pivastatin with_stat │       8 │  (31, 21) │     26 │     1 │  (2, 2) │     (4, 2) │   (1, 1) │    1.14021 │   0.0881271 │     0.00903408 │    nothing │     nothing │    nothing │   nothing │           26 │
│ Immune response to i │       9 │    (4, 1) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │  0.0115126 │ 0.000110838 │       3.772e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ Immune response to i │       7 │    (4, 1) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │ 0.00456527 │   5.9359e-5 │      2.0343e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ JAK-STAT 1 with_stat │      22 │  (19, 14) │     24 │     1 │  (2, 1) │    (2, 14) │   (1, 1) │   0.333888 │  0.00952475 │     0.00155373 │    nothing │     nothing │    nothing │   nothing │           24 │
│ QY with_states=false │      10 │   (13, 8) │    396 │     1 │  (4, 4) │    (10, 8) │  (11, 4) │   0.686769 │   0.0115483 │     0.00122763 │    nothing │     nothing │    nothing │   nothing │          342 │
│ HIV with_states=fals │      10 │    (8, 3) │     56 │     1 │  (3, 2) │     (3, 2) │   (1, 1) │   0.197885 │ 0.000247859 │       7.188e-5 │    nothing │     nothing │    nothing │   nothing │           28 │
│ LLW1987_io with_stat │       4 │    (3, 3) │     30 │     1 │  (2, 0) │     (2, 0) │   (2, 1) │   0.143245 │   3.0017e-5 │       1.454e-5 │    nothing │     nothing │    nothing │   nothing │           30 │
│ Ruminal lipolysis wi │       3 │    (3, 1) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │   0.141788 │   2.8912e-5 │      1.3022e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ Linear_compartment_h │       7 │    (5, 4) │   5412 │     1 │  (4, 4) │ (373, 363) │ (139, 1) │    336.287 │   0.0966693 │      0.0519126 │    nothing │     nothing │    nothing │   nothing │         5142 │
│ Pharm with_states=fa │       7 │  (27, 23) │     34 │     1 │  (1, 2) │     (1, 2) │   (1, 2) │    53.5016 │     0.84401 │        0.31273 │    nothing │     nothing │    nothing │   nothing │           20 │
│ SEIR2T with_states=f │       4 │    (7, 3) │     20 │     1 │  (1, 1) │     (1, 1) │   (1, 1) │ 0.00608748 │   8.6402e-5 │      3.0965e-5 │    nothing │     nothing │    nothing │   nothing │           18 │
│ SEUIR with_states=fa │       5 │   (14, 5) │     24 │     1 │  (2, 1) │     (2, 3) │   (1, 1) │  0.0154191 │ 0.000363824 │       8.813e-5 │    nothing │     nothing │    nothing │   nothing │           24 │
│ Immune response to i │       8 │    (4, 1) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │  0.0061465 │   6.0627e-5 │      2.2909e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ St with_states=false │       9 │  (18, 12) │    772 │     1 │  (4, 4) │    (11, 9) │  (26, 8) │    11.6043 │   0.0646193 │     0.00552757 │    nothing │     nothing │    nothing │   nothing │          662 │
│ Ovarian follicle pop │      22 │  (22, 18) │     34 │     1 │  (2, 1) │     (2, 9) │   (2, 1) │     306.34 │     1.84228 │       0.365265 │    nothing │     nothing │    nothing │   nothing │           34 │
│ Fujita with_states=f │      16 │  (15, 12) │     46 │     1 │  (1, 1) │     (4, 7) │   (3, 1) │    0.28889 │ 0.000847911 │     0.00018919 │    nothing │     nothing │    nothing │   nothing │           42 │
│ SIWR with extra outp │       7 │   (12, 9) │     20 │     1 │  (1, 1) │     (1, 1) │   (1, 1) │  0.0457402 │  0.00353208 │    0.000586929 │    nothing │     nothing │    nothing │   nothing │           18 │
│ Transfection_4State  │       5 │    (6, 3) │     20 │     1 │  (1, 1) │     (1, 1) │   (1, 1) │ 0.00493583 │   6.0907e-5 │      2.4542e-5 │    nothing │     nothing │    nothing │   nothing │           18 │
│ SIR 6 with_states=fa │       4 │    (2, 1) │     20 │     1 │  (1, 1) │     (1, 1) │   (1, 1) │ 0.00184782 │   1.8653e-5 │     1.21255e-5 │    nothing │     nothing │    nothing │   nothing │           18 │
│ HIV2_io with_states= │      10 │    (4, 4) │    372 │     1 │  (4, 2) │     (4, 2) │  (12, 1) │  0.0686153 │   7.9013e-5 │      3.9218e-5 │    nothing │     nothing │    nothing │   nothing │          278 │
│ Modified LV for test │       4 │    (2, 3) │     30 │     1 │  (2, 0) │     (2, 0) │   (2, 1) │ 0.00204316 │   2.0989e-5 │     1.18595e-5 │    nothing │     nothing │    nothing │   nothing │           30 │
│ CGV1990 with_states= │       9 │  (18, 10) │    132 │     1 │  (4, 3) │     (4, 3) │   (3, 2) │    0.26982 │  0.00303639 │    0.000603539 │    nothing │     nothing │    nothing │   nothing │           94 │
│ SEIR 36 ref with_sta │      13 │    (8, 5) │     20 │     1 │  (1, 0) │     (1, 3) │   (1, 1) │  0.0145667 │ 0.000506487 │      8.9661e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ Crauste_SI with_stat │      13 │    (4, 4) │     22 │     1 │  (1, 1) │     (2, 3) │   (1, 1) │  0.0119241 │ 0.000102357 │      4.0592e-5 │    nothing │     nothing │    nothing │   nothing │           18 │
│ SIWR original with_s │       7 │  (24, 11) │     24 │     1 │  (1, 2) │     (1, 2) │   (1, 1) │    2.15254 │    0.130931 │      0.0275474 │    nothing │     nothing │    nothing │   nothing │           20 │
│ Treatment_io with_st │       5 │    (7, 3) │     54 │     1 │  (2, 1) │     (4, 4) │   (3, 1) │ 0.00793072 │   6.6218e-5 │      2.6899e-5 │    nothing │     nothing │    nothing │   nothing │           54 │
│ cLV1 (2o) with_state │      15 │    (4, 4) │     26 │     1 │  (2, 2) │     (2, 2) │   (1, 1) │   0.238984 │ 0.000813357 │    0.000187775 │    nothing │     nothing │    nothing │   nothing │           26 │
│ generalizedLoktaVolt │       6 │    (3, 3) │     20 │     1 │  (1, 1) │     (1, 1) │   (1, 1) │ 0.00301457 │   5.0051e-5 │      2.4707e-5 │    nothing │     nothing │    nothing │   nothing │           18 │
│ Covid model (Gevertz │      10 │   (16, 8) │     26 │     1 │  (1, 0) │     (1, 3) │   (2, 1) │   0.124155 │  0.00811377 │     0.00112335 │    nothing │     nothing │    nothing │   nothing │           22 │
│ Immune response to i │       8 │    (3, 1) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │ 0.00737073 │   5.9895e-5 │      2.3292e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ Linear_compartment_h │       7 │    (5, 4) │   4388 │     1 │  (4, 2) │ (382, 373) │ (197, 1) │    200.292 │   0.0675165 │      0.0373318 │    nothing │     nothing │    nothing │   nothing │         4118 │
│ Immune response to i │       7 │    (3, 3) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │ 0.00285709 │   3.5137e-5 │      1.4706e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ Immune response to i │       6 │    (3, 3) │     18 │     1 │  (1, 0) │     (1, 0) │   (1, 1) │ 0.00242335 │   3.2648e-5 │      1.4627e-5 │    nothing │     nothing │    nothing │   nothing │           16 │
│ KD1999 with_states=f │      14 │   (10, 9) │     26 │     1 │  (2, 2) │     (2, 2) │   (1, 1) │    0.19712 │ 0.000356473 │      8.7097e-5 │    nothing │     nothing │    nothing │   nothing │           26 │
│ SLIQR with_states=fa │       6 │    (6, 4) │    716 │     1 │  (7, 6) │     (7, 6) │ (14, 11) │  0.0877847 │ 0.000116278 │      3.8949e-5 │    nothing │     nothing │    nothing │   nothing │          502 │
│ Chemical reaction ne │       6 │    (5, 4) │     34 │     1 │  (1, 2) │     (1, 2) │   (1, 2) │   0.012987 │ 0.000144237 │      5.0508e-5 │    nothing │     nothing │    nothing │   nothing │           20 │
│ MAPK model (5 output │      22 │   (10, 8) │     20 │     1 │  (1, 0) │    (1, 12) │   (1, 1) │    3.85126 │   0.0173041 │     0.00322229 │    nothing │     nothing │    nothing │   nothing │           16 │
│ SIR 24 with_states=f │       5 │   (11, 3) │     34 │     1 │  (2, 1) │     (2, 1) │   (2, 1) │  0.0140115 │ 0.000248774 │      8.0538e-5 │    nothing │     nothing │    nothing │   nothing │           34 │
│ Biohydrogenation_io  │       6 │    (6, 4) │     38 │     1 │  (2, 2) │     (2, 2) │   (2, 1) │ 0.00892917 │ 0.000114757 │        3.84e-5 │    nothing │     nothing │    nothing │   nothing │           38 │
└──────────────────────┴─────────┴───────────┴────────┴───────┴─────────┴────────────┴──────────┴────────────┴─────────────┴────────────────┴────────────┴─────────────┴────────────┴───────────┴──────────────┘

=#
