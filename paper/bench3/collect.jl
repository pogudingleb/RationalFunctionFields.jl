using StructuralIdentifiability, RationalFunctionFields, Nemo
import Base.Iterators
using PrettyTables

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

    table
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

    table
end

columns = [:name, :julia_funcs, :maple_funcs]
table1 = foo_julia_and_maple("table_funcs_julia_and_maple.md", columns)

columns = [:name, :julia_time, :maple_time]
table2 = foo_julia_and_maple("table_time_julia_and_maple.md", columns)

columns = [
        :name, :elems, :bytes, :max_deg,
        # :min_deg,
        :max_terms,
        # :min_terms,
        :vars,
        :min_deg_per_var,
        :min_gen_deg_per_var,
    ]
table3 = foo_gens_stats("table_input_stats.md", columns)

nothing
