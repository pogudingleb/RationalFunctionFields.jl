using Nemo

function param2()
R, (a,b,x0,x1,x2,x3,x4) = polynomial_ring(QQ, [:a,:b,:x0,:x1,:x2,:x3,:x4])

sys = [
    (b+1)*(x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0) - b,
   (a+b)*(2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1) - b,
   b*(x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2) - a,
   (a-1)*(2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3) - b,
   x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - a^3*big(2)^100
]
    sys, [x0,x1,x2,x3,x4], [a,b]
end
