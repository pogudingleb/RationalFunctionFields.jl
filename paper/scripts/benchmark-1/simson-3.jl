# Referenced in "Slimgb: Gröbner bases with slim polynomials", Michael Brickenstein, 2010. 
# Source code obtained from 
# https://github.com/symbolicdata/data/blob/a5cdf5366150bffba22a2c91a86a012dd23831af/XMLResources/IntPS/Geometry.Simson_3.xml
using Nemo

function simson3()
R, (u1, u2, u3, u4, x1, x2, x3, x4, x5, x6, x7, x8, x9) =
    polynomial_ring(QQ, [:u1, :u2, :u3, :u4, :x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9])

sys = [
    x1 * x8 - x2 * x8 - x1 * u2 + x2 * u1 - x9 * u1 + x9 * u2,
    -x2 * x6 + x7 * u2,
    -x1 * x4 + x5 * u1,
    x3^2 + u3^2 - 2 * u3 * u4,
    -x2 * x3 + x2 * x7 + x6 * u2 - u2 * u3,
    -x1 * x3 + x1 * x5 + x4 * u1 - u1 * u3,
    x1 * x3 - x2 * x3 - x1 * x9 + x2 * x9 - x8 * u1 + x8 * u2 + u1 * u3 - u2 * u3,
    x2^2 + u2^2 - 2 * u2 * u4,
    x1^2 + u1^2 - 2 * u1 * u4
]
    sys, [x1, x2, x3, x4, x5, x6, x7, x8, x9], [u1, u2, u3, u4]
end
