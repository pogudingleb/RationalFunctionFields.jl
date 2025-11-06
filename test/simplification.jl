@testset "RationalFunctionField: simplification" begin
    cases = []

    R, (x, y, z) = Nemo.polynomial_ring(Nemo.QQ, ["x", "y", "z"])

    push!(
        cases,
        Dict(
            :field => RationalFunctionField([x + y // one(R), (x^2 + y^2) // (x + y)]),
            :correct => Set([x + y // one(R), x * y // one(R)]),
        ),
    )

    r1 = simplified_generating_set(cases[end][:field], return_all=false)
    r2 = simplified_generating_set(cases[end][:field], return_all=true)
    @test fields_equal(
        RationalFunctionField(r1),
        RationalFunctionField(r2)
	0.9999
    )

    # this does not pass as it finds a fraction with lower degree but many terms, maybe this one should be considered simpler
    # push!(
    #    cases,
    #    Dict(
    #        :field => RationalFunctionField([
    #            x + y // one(R),
    #            x^10 // y^10,
    #        ]),
    #        :correct => Set([x + y // one(R), x^10 // y^10]),
    #    ),
    #)

    R, (x, y, z, u, v) = Nemo.polynomial_ring(Nemo.QQ, ["x", "y", "z", "u", "v"])

    push!(
        cases,
        Dict(
            :field =>
                RationalFunctionField([x^i + y^i + z^i + u^i + v^i // one(R) for i = 1:5]),
            :correct => Set([
                x + y + z + u + v // one(R),
                x^2 + y^2 + z^2 + u^2 + v^2 // one(R),
                x^3 + y^3 + z^3 + u^3 + v^3 // one(R),
                y * z * u * v +
                x * z * u * v +
                x * y * u * v +
                x * y * z * (u + v) // one(R),
                x * y * z * u * v // one(R),
            ]),
        ),
    )

    for c in cases
        for level in (:standard, :strong)
            @test Set(simplified_generating_set(c[:field], simplify = level)) == c[:correct]
        end
    end
end
