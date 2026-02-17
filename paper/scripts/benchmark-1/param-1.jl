using Nemo

function param1()
    R, (a, b, x0, x1, x2, x3, x4) = polynomial_ring(
        QQ,
        [:a, :b, :x0, :x1, :x2, :x3, :x4],
        internal_ordering = :degrevlex,
    )
    sys = [
        x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0 - b,
        2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1,
        x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2,
        2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3,
        b*x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - a,
    ]
    sys, [x0, x1, x2, x3, x4], [a, b]
end
