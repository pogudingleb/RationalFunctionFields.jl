@testset "MQS raw ideal generators" begin
    R, (a, b, c) = QQ["a", "b", "c"]
    funcs = [[R(1), R(1)]]
    mqs = IdealMQS(funcs)
    sys, _ = fractionfree_generators_raw(mqs)
    @test sys == [parent(sys[1])(0)]
    
    funcs = [[R(1), a], [b, a + c, c^2, b]]
    mqs = IdealMQS(funcs)
    sys, indets, params = fractionfree_generators_raw(mqs)
    (t1,y1,y2,y3) = indets
    (a,b,c) = params
    @test sys == [
        y1 - a,
        -(a+c)*y2 + b*(y1+y3),
        -c^2*y2 + b*y3^2,
        parent(t1)(0),
        t1*y2 - 1
    ]
end
