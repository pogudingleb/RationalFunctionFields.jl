import Pkg; Pkg.activate(joinpath(@__DIR__, "..", "env"))
using StructuralIdentifiability, RationalFunctionFields, Nemo
import Base.Iterators
using PrettyTables
using Latexify, BenchmarkTools

function myprettymemory(b)
    if b < 1024
        return string(b, "~bytes")
    elseif b < 1024^2
        value, units = round(b / 1024.0, digits=1), "KB"
    elseif b < 1024^3
        value, units = round(b / 1024.0^2, digits=1), "MB"
    else
        value, units = round(b / 1024.0^3, digits=1), "GB"
    end
    return string(value, "~", units)
end

# get `benchmarks`
if !isdefined(Main, :benchmarks)
include(joinpath(dirname(dirname(pathof(StructuralIdentifiability))), "benchmarking", "benchmarks.jl"))
end

skipped = ["NFkB", "Covid2", "Akt", "TumorPillis", "TumorHu", "LeukaemiaLeon2021", "MAPK_5out_bis", "cLV2", "QWWC"]
prefix = "results"

function get_julia_result(name)
    output_filepath = joinpath(@__DIR__, prefix, "simplify", string(name, "_output.txt"))
    result = Dict{Any,Any}(
        :julia_funcs => nothing, 
        :julia_time => nothing,
        :julia_max_deg => nothing,
        :julia_min_deg => nothing,
        :julia_max_terms => nothing,
        :julia_min_terms => nothing
    )    
    !isfile(output_filepath) && return result
    open(output_filepath, "r") do io
        last_line = nothing
        for line in Iterators.reverse(eachline(io))
            if startswith(line, ">>>")
                result[:julia_time] = chopprefix(line, ">>>")
            end
            if startswith(line, "simple funcs")
                result[:julia_funcs] = chopprefix(last_line, "AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}")
            end
            for (str, key) in zip(["max deg = ", "min deg = ", "max terms = ", "min terms = "], [:julia_max_deg, :julia_min_deg, :julia_max_terms, :julia_min_terms])
                if startswith(line, str)
                    result[key] = chopprefix(line, str)
                end
            end
            last_line = line
            !any(nothing .== collect(values(result))) && break
        end
    end
    return result
end

function get_maple_input_gens(name)
    output_filepath = joinpath(@__DIR__, prefix, "maple_simplify", string(name, "_gens.mpl"))
    result = Dict{Any,Any}(
        :maple_input_gens => nothing, 
    )
    !isfile(output_filepath) && return result
    open(output_filepath, "r") do io
        content = read(io, String)
        content = strip(content)
        content = chopprefix(content, "gens := {")
        content = chopsuffix(content, "}:")
        result[:maple_input_gens] = content
    end
    return result
end

function get_maple_result(name)
    output_filepath = joinpath(@__DIR__, prefix, "maple_simplify", string(name, "_run_simplify_output.txt"))
    result = Dict{Any,Any}(
        :maple_funcs => nothing, 
        :maple_time => nothing,
    )
    !isfile(output_filepath) && return result
    open(output_filepath, "r") do io
        last_line = nothing
        for line in Iterators.reverse(eachline(io))
            if occursin("time=", line) && result[:maple_time] == nothing
                result[:maple_time] = split(line, "time=")[2]
            end
            if startswith(line, "\">>>\"")
                result[:maple_funcs] = chopprefix(line, "\">>>\", ")
            end
            !any(nothing .== collect(values(result))) && break
        end
    end
    if strip(string(result[:maple_funcs])) == "{gens}"
        result[:maple_funcs] = nothing
    end
    return result
end

function get_gens_stats(name)
    output_filepath = joinpath(@__DIR__, prefix, "input_stats", string(name, "_gens_stats.txt"))
    result = Dict{Any,Any}(
        :max_deg => nothing,
        :min_deg => nothing,
        :max_terms => nothing,
        :min_terms => nothing,
        :elems => nothing,
        :bytes => nothing,
        :bytes_dennums => nothing,
        :min_gen_deg_per_var => nothing,
        :min_deg_per_var => nothing,
        :vars => nothing
    )    
    !isfile(output_filepath) && return result
    open(output_filepath, "r") do io
        last_line = ""
        for line in Iterators.reverse(eachline(io))
            for (str, key) in zip(["max deg(...) = ", "min deg(...) = ", "max terms(...) = ", "min terms(...) = ", "length(...) = ", "min gen deg per var = ", "min deg per var = ", "vars = "], [:max_deg, :min_deg, :max_terms, :min_terms, :elems, :min_gen_deg_per_var, :min_deg_per_var, :vars])
                if startswith(line, str)
                    result[key] = chopprefix(line, str)
                end
            end
            for (str, key) in zip(["format [num/den]", "format [[den, nums]]"], [:bytes, :bytes_dennums])
                if startswith(line, str)
                    result[key] = chopprefix(last_line, "length(string(...)) = ")
                end
            end
            last_line = line
            !any(nothing .== collect(values(result))) && break
        end
    end
    return result
end

function get_independence_result(name)
    output_filepath = joinpath(@__DIR__, prefix, "independence", string(name, "_output.txt"))
    result = Dict{Any,Any}(
        :independent => nothing,
    )    
    !isfile(output_filepath) && return result
    open(output_filepath, "r") do io
        last_line = ""
        for line in Iterators.reverse(eachline(io))
            for (str, key) in zip(["independent: "], [:independent])
                if startswith(line, str)
                    result[key] = chopprefix(line, str)
                end
            end
            last_line = line
            !any(nothing .== collect(values(result))) && break
        end
    end
    return result
end

function foo_get_maple_input_gens()
    results = []
    for name in keys(benchmarks)
        if any(s -> occursin(s, string(name)), skipped)
            @error "Skipping $name"
            continue
        end
        maple_result = get_maple_input_gens(name)
        common = Dict(:name => name)
        result = merge(common, maple_result)
        push!(results, result)
    end

    results = sort(results, by=x -> x[:name])
    results
end

function foo_julia_and_maple(filename, columns)
    results = []
    for name in keys(benchmarks)
        if any(s -> occursin(s, string(name)), skipped)
            @error "Skipping $name"
            continue
        end
        julia_result = get_julia_result(name)
        maple_result = get_maple_result(name)
        common = Dict(:name => name)
        result = merge(common, julia_result, maple_result)
        push!(results, result)
    end

    results = sort(results, by=x -> x[:name])

    table = Matrix{Any}(undef, length(results), length(columns))
    for (i, result) in enumerate(results)
        table[i, :] = [result[c] for c in columns]
    end

    result_path = joinpath(@__DIR__, filename)
    open(result_path, "w") do io
        pretty_table(
            io,
            table, 
            column_labels = columns,
            backend = :markdown
        )
        println("Result written to $result_path")
    end

    results
end

function foo_julia_original_generators(columns)
    results = []
    for name in keys(benchmarks)
        if any(s -> occursin(s, string(name)), skipped)
            @error "Skipping $name"
            continue
        end
        common = Dict(:name => name)
        initial_funcs = nothing
        output_filepath = joinpath(@__DIR__, prefix, "input_stats", string(name, "_initial_funcs.txt"))
        if ispath(output_filepath)
            if filesize(output_filepath) < 2^21 # 2 MB
                open(output_filepath, "r") do io
                    content = read(io, String)
                    content = chopprefix(content, "AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}")
                    initial_funcs = content
                end
            end
        end
        result = merge(common, Dict(:original_funcs => initial_funcs))
        push!(results, result)
    end

    results = sort(results, by=x -> x[:name])
    
    results
end

function foo_gens_stats(filename, columns)
    results = []
    for name in keys(benchmarks)
        if any(s -> occursin(s, string(name)), skipped)
            @error "Skipping $name"
            continue
        end
        gens_stats_result = get_gens_stats(name)
        common = Dict(:name => name)
        result = merge(common, gens_stats_result)
        push!(results, result)
    end

    results = sort(results, by=x -> x[:name])
    
    table = Matrix{Any}(undef, length(results), length(columns))
    for (i, result) in enumerate(results)
        table[i, :] = [result[c] for c in columns]
    end

    result_path = joinpath(@__DIR__, filename)
    open(result_path, "w") do io
        pretty_table(
            io,
            table, 
            column_labels = columns,
            backend = :markdown
        )
        println("Result written to $result_path")
    end

    results
end

function foo_independence(filename, columns)
    results = []
    for name in keys(benchmarks)
        if any(s -> occursin(s, string(name)), skipped)
            @error "Skipping $name"
            continue
        end
        independence_result = get_independence_result(name)
        common = Dict(:name => name)
        result = merge(common, independence_result)
        push!(results, result)
    end

    results = sort(results, by=x -> x[:name])
    
    table = Matrix{Any}(undef, length(results), length(columns))
    for (i, result) in enumerate(results)
        table[i, :] = [result[c] for c in columns]
    end

    result_path = joinpath(@__DIR__, filename)
    open(result_path, "w") do io
        pretty_table(
            io,
            table, 
            column_labels = columns,
            backend = :markdown
        )
        println("Result written to $result_path")
    end

    results
end

function paginate(words; chars_per_page = 40, sep = ",")
    pages = []
    current_page = ""
    for word in words
        if length(current_page) + length(word) + length(sep) > chars_per_page
            push!(pages, current_page)
            current_page = ""
        end
        current_page *= word * sep
    end
    push!(pages, current_page)
    pages
end

function nicify(result, f; do_paginate=true)
    f = string(f)
    f = strip(f)
    f = chopprefix(f, "[")
    f = chopsuffix(f, "]")
    f = replace(f, r"__" => "_")
    f = replace(f, "gamma" => "gama")
    # f = replace(f, "//" => "/")
    f = split(f, ",")
    f = map(strip, f)
    f = sort(f, by=c -> (length(c), c))
    f = map(ff -> replace(ff, r"reaction\_([0-9])\_k([0-9])" => s"r_\1_\2"), f)
    f = map(strip, f)
    if do_paginate
        f = paginate(f, chars_per_page = 80, sep = ", ")
    end
    f = map(strip, f)
    f = filter(!isempty, f)
    f = map(x -> strip(x, ','), f)
    println("before latexify: ", f)
    f = map(x -> join(String.(latexify.(split(x, ", "), mult_symbol="")), ", "), f)
    println("after String(latex()): ", f)
    if do_paginate
        f = join(f, ",\\\\")
    else
        f = join(f, ", ")
    end
    f = replace(f, r"\\frac" => "\\dfrac")
    f = replace(f, r"Mar" => "M")
    f = replace(f, r"Ninv" => "N_{\\operatorname{inv}}")
    f = replace(f, r"siga([0-9]+)" => s"s_{\1}")
    f = replace(f, "gama" => "gamma")
    f = replace(f, "alpa" => "alpha")
    f = replace(f, r"([a-zA-Z]+)([0-9]+)" => s"\1_{\2}")
    greek = ["mu", "rho", "beta", "gamma", "alpha", "delta", "phi", "tau"]
    for x in vcat(greek, uppercasefirst.(greek))
        f = replace(f, Regex("([\$ _{}])$(x)([\$ _{}])") => SubstitutionString("\\1\\\\$(lowercase(x))\\2"))
    end
    f = replace(f, r"_{([0-9])\\_([0-9])}" => s"_{\1\2}")
    for _ in 1:10 f = replace(f, "  " => " ") end
    f = replace(f, ", " => ",")
    f = replace(f, "\$ " => "\$")
    f = replace(f, "\$- " => "\$-")
    f = replace(f, " - " => "-")
    f = replace(f, " + " => "+")
    if string(result[:name]) in ["Bruno2016", "Covid3", "EAIHRD", "HDNL", "HIV", "HIV2", "Influenza_MB1", "Influenza_MB3", "Influenza_MD1", "Influenza_MD2", "Influenza_MD3", "Influenza_MD4", "SEAIJRC", "SIR19", "SIR21", "SLIQR", "Transfection", "Treatment"]
        f = replace(f, " " => "\\mspace{2mu} ")
    end
    f = replace(f, "," => ", ")
    f = String(f)
    if string(result[:name]) in ["Lincomp1", "Lincomp2"]
        newf = []
        f = split(f, ",")
        for (i, el) in enumerate(f)
            trms = reduce(vcat, map(e -> split(e, "-"), split(el, "+")))
            maxtrms = 5
            if length(el) > 200 # length(trms) > maxtrms
                push!(newf, string(join(trms[1:maxtrms], "+"), " + \\text{$(length(trms) - maxtrms) more terms}\$"))
            else
                push!(newf, el)
            end
        end
        newf = join(newf, ",")
        f = newf
    end
    f
end

function make_nice_latex_table(filename, columns, results)

    for result in results
        println(result[:name])
        f = result[:julia_funcs]
        f = nicify(result, f)
        rowspace = "\\setlength{\\extrarowheight}{$(occursin("frac", f) ? 8 : 0)pt}"
        f = "{\\footnotesize$(rowspace)\\begin{tabular}{l}$(f).\\end{tabular}}"
        result[:julia_funcs] = f

        f = result[:maple_input_gens]
        if !isnothing(f) && length(f) < 500
            f = nicify(result, f)
            rowspace = "\\setlength{\\extrarowheight}{$(occursin("frac", f) ? 8 : 0)pt}"
            f = "{\\footnotesize$(rowspace)\\begin{tabular}{l}$(f).\\end{tabular}}"
            result[:maple_input_gens] = f
        else
            result[:maple_input_gens] = nothing
        end

        f = result[:name]
        f = string(f)
        f = replace(f, "_" => "\\_")
        result[:name] = f

        f = result[:vars]
        if !isnothing(f)
            f = nicify(result, f; do_paginate=false)
            result[:vars] = f
        end
    end

    results = filter(result -> !(result[:name] in ["QY", "St"]), results)

    table = Matrix{Any}(undef, length(results), length(columns))
    for (i, result) in enumerate(results)
        table[i, :] = [result[c] for c in columns]
    end

    columns = map(x -> replace(string(x), "_" => "\\_"), columns)

    result_path = joinpath(@__DIR__, filename)
    open(result_path, "w") do io
        println(io, """
        \\begin{enumerate}
            """)
        for result in results
            name = result[:name]
            julia_funcs = result[:julia_funcs]
            input_funcs = get(result, :maple_input_gens, nothing)
            independent = get(result, :independent, nothing)
            elems = get(result, :elems, nothing)
            vars = get(result, :vars, nothing)
            bytes = get(result, :bytes, nothing)
            if isnothing(bytes) bytes = "0" end
            max_deg = get(result, :max_deg, nothing)
            println(io, "\\item \\myexample{$(name)} \\cite{}.\n")
            println(io, "Original generating set information: $(length(split(vars, ","))) indeterminates; $elems non-constant functions; maximal total degrees of numerator and denominator are~\$$max_deg\$; $(myprettymemory(parse(Int, bytes))) total in string representation.\n")
            println(io, "Indeterminates: $vars.\n")
            if input_funcs !== nothing
            println(io, "Original generating set: \n\n\$$(input_funcs)\$\n")
            else
            println(io, "Original generating set: too large to be listed.\n")
            end
            println(io, "Result of our algorithm: \n\n\$$(julia_funcs)\$\n")
            println(io, "Result is algebraically independent over \$\\mathbb{C}\$: $(independent in ["true", nothing] ? "yes" : "no").")
            println(io)
        end
        println(io, """
        \\end{enumerate}
        """)
        println("Result written to $result_path")
    end

    table
end

function merge_by_name(results_lists...)
    results = []
    for result in results_lists[1]
        name = result[:name]
        for ress in results_lists[2:end]
            found = findall(res -> res[:name] == name, ress)
            isempty(found) && continue
            result = merge(result, only(ress[found]))
        end
        push!(results, result)
    end
    results
end

function print_huge_data_to_data_1(
        results; 
        filename = nothing,
        column = nothing,
        size_limit = 2^19   # half a MB
    )
    path = joinpath(dirname(@__DIR__), "data-1")
    @assert ispath(path)
    for result in results
        model_path = joinpath(path, replace(string(result[:name]), " " => "_", "." => "_"))
        !ispath(model_path) && mkdir(model_path)
        content = result[column]
        if content == nothing || length(content) > size_limit
            continue
        end
        open(joinpath(model_path, filename), "w") do io
            println(io, content)
        end
    end
end

function main()
    task = :all
    if !isempty(ARGS)
        task = ARGS[1]
    end
    
    columns = [:name, :julia_funcs]
    results1 = foo_julia_and_maple("table_funcs_julia_and_maple_fan.md", columns)
    print_huge_data_to_data_1(results1, column=:julia_funcs, filename="our_output.txt")
    
    columns = [:name, :original_funcs]
    results11 = foo_julia_original_generators(columns)
    print_huge_data_to_data_1(results11, column=:original_funcs, filename="original_generators.txt")

    rat = filter(res -> res[:original_funcs] == nothing || occursin("/", res[:original_funcs]), results11);
    poly_to_rat = filter(res -> !occursin("/", results1[findfirst(res1 -> res1[:name] == res[:name], results1)][:julia_funcs]), filter(res -> res[:original_funcs] != nothing, rat))
    @error "" length(poly_to_rat)
    
    columns = [:name, :julia_time, :maple_time]
    table2 = foo_julia_and_maple("table_time_julia_and_maple_fan.md", columns)
    
    columns = [
            :name, :elems, :bytes, :max_deg,
            # :min_deg,
            :max_terms,
            # :min_terms,
            :vars,
            :min_deg_per_var,
            :min_gen_deg_per_var,
        ]
    results2 = foo_gens_stats("table_input_stats.md", columns)
    
    results3 = foo_independence("table_independence.md", [:name, :independent])
    
    results4 = foo_get_maple_input_gens()
    
    results = merge_by_name(results1,results2,results3,results4)
    table11 = make_nice_latex_table("table.tex", columns, results)

end


main()
