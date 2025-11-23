@testset "MQS raw ideal generators" begin
    # name collision
    _, (x, y1) = QQ["x", "y1"]
    funcs = [[x, y1]]
    mqs = IdealMQS(funcs)
    @test_throws AssertionError fractionfree_generators_raw(mqs)

    R, (a, b, c, x) = QQ["a", "b", "c", "x(t)"]
    funcs = [[R(1), R(1)]]
    mqs = IdealMQS(funcs)
    sys, _ = fractionfree_generators_raw(mqs)
    @test sys == [parent(sys[1])(0)]

    funcs = [[R(1), x], [b, a + c, c^2, b]]
    mqs = IdealMQS(funcs)
    sys, indets, params = fractionfree_generators_raw(mqs)
    (t1,y1,y2,y3,y4) = indets
    (a,b,c,x) = params
    @test sys == [
        y4 - x,
        -(a+c)*y2 + b*(y1+y3),
        -c^2*y2 + b*y3^2,
        parent(t1)(0),
        t1*y2 - 1
    ]

    K = Nemo.Native.GF(2^30+3)
    RationalFunctionFields.reduce_mod_p!(mqs, K)
    p = rand(K, length(params))
    mqs_spec_mod_p = RationalFunctionFields.specialize_mod_p(mqs, p)
    R = parent(mqs_spec_mod_p[1])
    poly_mod_p(f, K) = map_coefficients(c -> K(numerator(c)) // K(denominator(c)), f)
    gens_mod_p = map(f -> poly_mod_p(f, K), sys)
    params_mod_p = map(f -> poly_mod_p(f, K), params)
    indets_mod_p = map(f -> poly_mod_p(f, K), indets)
    gens_spec_mod_p = map(f -> evaluate(f, params_mod_p, p), gens_mod_p)
    gens_spec_mod_p = map(f -> RationalFunctionFields.parent_ring_change(f, R, matching=:byname), gens_spec_mod_p)
    @test mqs_spec_mod_p == gens_spec_mod_p
end
