using Nemo
using RationalFunctionFields
using RationalFunctionFields: field_contains_mod_p

R, (a1, a2, a3, b1, b2, b3) = polynomial_ring(QQ, ["a1", "a2", "a3", "b1", "b2", "b3"])


gens = [
(a2-a3)*(a1-a3)*(a1-a2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(a1*b2^2-a1*b3^2-a2*b1^2+a2*b3^2+a3*b1^2-a3*b2^2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(a1*b1*b2-a1*b1*b3-a2*b1*b2+a2*b2*b3+a3*b1*b3-a3*b2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(b1^2*b2+b1^2*b3+b1*b2^2-6*b1*b2*b3+b1*b3^2+b2^2*b3+b2*b3^2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(2*b1^3-b1^2*b2-b1^2*b3-b1*b2^2-b1*b3^2+2*b2^3-b2^2*b3-b2*b3^2+2*b3^3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(a1*a2*b1+a1*a2*b2-2*a1*a2*b3+a1*a3*b1-2*a1*a3*b2+a1*a3*b3-2*a2*a3*b1+a2*a3*b2+a2*a3*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(2*a1^2*b1-a1^2*b2-a1^2*b3-a2^2*b1+2*a2^2*b2-a2^2*b3-a3^2*b1-a3^2*b2+2*a3^2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(a1*b1^2*b2+a1*b1^2*b3-2*a1*b1*b2*b3+a2*b1*b2^2-2*a2*b1*b2*b3+a2*b2^2*b3-2*a3*b1*b2*b3+a3*b1*b3^2+a3*b2*b3^2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(a1*b1*b2-a1*b1*b3-a1*b2^2+a1*b3^2+a2*b1^2-a2*b1*b2+a2*b2*b3-a2*b3^2-a3*b1^2+a3*b1*b3+a3*b2^2-a3*b2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(a1*b1*b2^2+a1*b1*b3^2-a1*b2^3-a1*b3^3-a2*b1^3+a2*b1^2*b2+a2*b2*b3^2-a2*b3^3-a3*b1^3+a3*b1^2*b3-a3*b2^3+a3*b2^2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(a1+a2+a3)*(a1*b1*b2-a1*b1*b3-a2*b1*b2+a2*b2*b3+a3*b1*b3-a3*b2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(2*a1*b1^3-a1*b1^2*b2-a1*b1^2*b3+a1*b2^3-a1*b2^2*b3-a1*b2*b3^2+a1*b3^3+a2*b1^3-a2*b1^2*b3-a2*b1*b2^2-a2*b1*b3^2+2*a2*b2^3-a2*b2^2*b3+a2*b3^3+a3*b1^3-a3*b1^2*b2-a3*b1*b2^2-a3*b1*b3^2+a3*b2^3-a3*b2*b3^2+2*a3*b3^3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
(2*a1^2*b1^2+a1^2*b1*b2+a1^2*b1*b3-a1*a2*b1^2-a1*a2*b1*b3-a1*a2*b2^2-a1*a2*b2*b3-a1*a3*b1^2-a1*a3*b1*b2-a1*a3*b2*b3-a1*a3*b3^2+a2^2*b1*b2+2*a2^2*b2^2+a2^2*b2*b3-a2*a3*b1*b2-a2*a3*b1*b3-a2*a3*b2^2-a2*a3*b3^2+a3^2*b1*b3+a3^2*b2*b3+2*a3^2*b3^2)//(b3+b1+b2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(a2-a3)*(a1-a3)*(a1-a2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(a1*b2^2-a1*b3^2-a2*b1^2+a2*b3^2+a3*b1^2-a3*b2^2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-2*(a1*b1*b2-a1*b1*b3-a2*b1*b2+a2*b2*b3+a3*b1*b3-a3*b2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(a1*b1*b2-a1*b1*b3-a2*b1*b2+a2*b2*b3+a3*b1*b3-a3*b2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
2*(a1*b1*b2-a1*b1*b3-a2*b1*b2+a2*b2*b3+a3*b1*b3-a3*b2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(b1^2*b2+b1^2*b3+b1*b2^2-6*b1*b2*b3+b1*b3^2+b2^2*b3+b2*b3^2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(2*b1^3-b1^2*b2-b1^2*b3-b1*b2^2-b1*b3^2+2*b2^3-b2^2*b3-b2*b3^2+2*b3^3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(a1^2*b1*b2+a1^2*b1*b3-2*a1*a2*b1*b2-2*a1*a3*b1*b3+a2^2*b1*b2+a2^2*b2*b3-2*a2*a3*b2*b3+a3^2*b1*b3+a3^2*b2*b3)//(b3+b1+b2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(a1*b1*b2-a1*b1*b3-a1*b2^2+a1*b3^2+a2*b1^2-a2*b1*b2+a2*b2*b3-a2*b3^2-a3*b1^2+a3*b1*b3+a3*b2^2-a3*b2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-1//2*(a1*b2+a1*b3+a2*b1+a2*b3+a3*b1+a3*b2)*(a1*b2-a1*b3-a2*b1+a2*b3+a3*b1-a3*b2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
1//2*(a1*b2+a1*b3+a2*b1+a2*b3+a3*b1+a3*b2)*(a1*b2-a1*b3-a2*b1+a2*b3+a3*b1-a3*b2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(a1+a2+a3)*(a1*b1*b2-a1*b1*b3-a2*b1*b2+a2*b2*b3+a3*b1*b3-a3*b2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(a1^2*a2*b2-a1^2*a2*b3-a1^2*a3*b2+a1^2*a3*b3+a1*a2^2*b1-a1*a2^2*b3+a1*a3^2*b1-a1*a3^2*b2-a2^2*a3*b1+a2^2*a3*b3-a2*a3^2*b1+a2*a3^2*b2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(a1^2*a2*b1-a1^2*a2*b2+a1^2*a3*b1-a1^2*a3*b3-a1*a2^2*b1+a1*a2^2*b2-a1*a3^2*b1+a1*a3^2*b3+a2^2*a3*b2-a2^2*a3*b3-a2*a3^2*b2+a2*a3^2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-(a1*b1*b2^2-4*a1*b1*b2*b3+a1*b1*b3^2+a1*b2^2*b3+a1*b2*b3^2+a2*b1^2*b2+a2*b1^2*b3-4*a2*b1*b2*b3+a2*b1*b3^2+a2*b2*b3^2+a3*b1^2*b2+a3*b1^2*b3+a3*b1*b2^2-4*a3*b1*b2*b3+a3*b2^2*b3)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-1//2*(a1^2*b2^2-a1^2*b3^2-2*a1*a2*b1^2+2*a1*a2*b1*b3+2*a1*a2*b2^2-2*a1*a2*b2*b3+2*a1*a3*b1^2-2*a1*a3*b1*b2+2*a1*a3*b2*b3-2*a1*a3*b3^2-a2^2*b1^2+a2^2*b3^2+2*a2*a3*b1*b2-2*a2*a3*b1*b3-2*a2*a3*b2^2+2*a2*a3*b3^2+a3^2*b1^2-a3^2*b2^2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
1//2*(a1^2*b2^2-a1^2*b3^2-2*a1*a2*b1^2+2*a1*a2*b1*b3+2*a1*a2*b2^2-2*a1*a2*b2*b3+2*a1*a3*b1^2-2*a1*a3*b1*b2+2*a1*a3*b2*b3-2*a1*a3*b3^2-a2^2*b1^2+a2^2*b3^2+2*a2*a3*b1*b2-2*a2*a3*b1*b3-2*a2*a3*b2^2+2*a2*a3*b3^2+a3^2*b1^2-a3^2*b2^2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-1//2*(a1^2*b2^2+2*a1^2*b2*b3+a1^2*b3^2+2*a1*a2*b1*b2-2*a1*a2*b1*b3-2*a1*a2*b2*b3-2*a1*a2*b3^2-2*a1*a3*b1*b2+2*a1*a3*b1*b3-2*a1*a3*b2^2-2*a1*a3*b2*b3+a2^2*b1^2+2*a2^2*b1*b3+a2^2*b3^2-2*a2*a3*b1^2-2*a2*a3*b1*b2-2*a2*a3*b1*b3+2*a2*a3*b2*b3+a3^2*b1^2+2*a3^2*b1*b2+a3^2*b2^2)//(b3+b1+b2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
-1//2*(a1^2*b2^2+2*a1^2*b2*b3+a1^2*b3^2-2*a1*a2*b1^2-2*a1*a2*b1*b2-2*a1*a2*b2^2+2*a1*a2*b3^2-2*a1*a3*b1^2-2*a1*a3*b1*b3+2*a1*a3*b2^2-2*a1*a3*b3^2+a2^2*b1^2+2*a2^2*b1*b3+a2^2*b3^2+2*a2*a3*b1^2-2*a2*a3*b2^2-2*a2*a3*b2*b3-2*a2*a3*b3^2+a3^2*b1^2+2*a3^2*b1*b2+a3^2*b2^2)//(b3+b1+b2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
1//2*(2*a1^2*b1*b2+2*a1^2*b1*b3+a1^2*b2^2+2*a1^2*b2*b3+a1^2*b3^2-2*a1*a2*b1*b2-2*a1*a2*b1*b3-2*a1*a2*b2*b3-2*a1*a2*b3^2-2*a1*a3*b1*b2-2*a1*a3*b1*b3-2*a1*a3*b2^2-2*a1*a3*b2*b3+a2^2*b1^2+2*a2^2*b1*b2+2*a2^2*b1*b3+2*a2^2*b2*b3+a2^2*b3^2-2*a2*a3*b1^2-2*a2*a3*b1*b2-2*a2*a3*b1*b3-2*a2*a3*b2*b3+a3^2*b1^2+2*a3^2*b1*b2+2*a3^2*b1*b3+a3^2*b2^2+2*a3^2*b2*b3)//(b3+b1+b2)//(2*a1*b1-a1*b2-a1*b3-a2*b1+2*a2*b2-a2*b3-a3*b1-a3*b2+2*a3*b3),
a1+a2+a3
]

function td(f)
    return total_degree(numerator(f // one(R))) + total_degree(denominator(f // one(R)))
end

println("There are $(length(gens)) original generators, their total degrees are $(map(td, gens))")
sg = simplified_generating_set(RationalFunctionField(gens), enforce_minimality = true)
println(sg)
