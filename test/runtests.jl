using RationalFunctionFields

using Test
using TestSetExtensions

using Nemo

function get_test_files()
    result = Vector{String}()
    for (dir, _, files) in walkdir("./")
        for fname in files
            if fname != "runtests.jl" && endswith(fname, ".jl")
                push!(result, dir * "/" * fname)
            end
        end
    end
    return result
end

@info "Testing started"

all_tests = get_test_files()
if !isempty(ARGS)
    all_tests = ARGS
end

@time @testset "All the tests" verbose = true begin
    for test_file in all_tests
            include(test_file)
    end
end
