using RationalFunctionFields
using Nemo

R, (a10, a01, a20, a11, a02, b10, b01, b20, b11, b02) = polynomial_ring(QQ, ["a10", "a01", "a20", "a11", "a02", "b10", "b01", "b20", "b11", "b02"])

gens = [
(a01-b01)*(a01+b01)//(a01*a02*b10-a01*a11*b01+a01*b01*b11-a10*b01*b02),
(a01^2*a02-b01^2*b02)//(a02-b02)//(a02+b02),
(a02*b11-a11*b02)//(a02-b02),
(a02^3*a20-a02^2*a11^2-b02^3*b20+b02^2*b11^2)//(a02-b02)//(a02+b02),
(a01^2*b11-a01*a10*b02+a02*b01*b10-a11*b01^2)//(a01*a02*b10-a01*a11*b01+a01*b01*b11-a10*b01*b02),
(a01^2*b02*b11-a01*a10*b02^2+a02^2*b01*b10-a02*a11*b01^2)//(a02-b02)//(a02+b02),
(a02^3*b20-2*a02^2*a11*b11+a02*a11^2*b02-a02*b02*b11^2+2*a11*b02^2*b11-a20*b02^3)//(a02-b02)//(a02+b02),
(a01^2*a02*b11^2-2*a01*a02*a10*b02*b11-a02^2*b02*b10^2+a02*a10^2*b02^2+2*a02*a11*b01*b02*b10-a11^2*b01^2*b02)//(a02-b02)//(a02+b02),
(a02^2+b02^2)//a02//b02,
(a01*a02^2*b01*b11+a01*a02*b02^2*b10-a01*a11*b01*b02^2-a02^2*a10*b01*b02)//a02//b02//(a01*a02*b10-a01*a11*b01+a01*b01*b11-a10*b01*b02),
(a01*a02^3*b10-a01*a02^2*a11*b01+a01*b01*b02^2*b11-a10*b01*b02^3)//a02//b02//(a01*a02*b10-a01*a11*b01+a01*b01*b11-a10*b01*b02),
-(a01^2*b02-a02*b01^2)//(a02-b02)//(a02+b02),
-(a02*b11-a11*b02)//(a02-b02),
-(a02^2*b02*b20-a02^2*b11^2-a02*a20*b02^2+a11^2*b02^2)//(a02-b02)//(a02+b02),
-(a01^2*a11-a01*a02*a10-b01^2*b11+b01*b02*b10)//(a01*a02*b10-a01*a11*b01+a01*b01*b11-a10*b01*b02),
-(a01^2*a02*b11-a01*a02*a10*b02+a02*b01*b02*b10-a11*b01^2*b02)//(a02-b02)//(a02+b02),
-(a01^2*b02*b11^2-2*a01*a10*b02^2*b11-a02^3*b10^2+2*a02^2*a11*b01*b10-a02*a11^2*b01^2+a10^2*b02^3)//(a02-b02)//(a02+b02),
-(a01^2*a11*b11-a01*a02*a10*b11-a01*a10*a11*b02+a02*a10^2*b02+a02*b01*b10*b11-a02*b02*b10^2-a11*b01^2*b11+a11*b01*b02*b10)//(a01*a02*b10-a01*a11*b01+a01*b01*b11-a10*b01*b02),
-a02*b02*(a02*a20-a11^2-b02*b20+b11^2)//(a02-b02)//(a02+b02),
-(a02-b02)*(a02+b02)*(a02*b10-a11*b01)*(a01*b11-a10*b02)//a02//b02//(a01*a02*b10-a01*a11*b01+a01*b01*b11-a10*b01*b02),
-a01*b01*(a02-b02)*(a02+b02)//a02//b02//(a01*a02*b10-a01*a11*b01+a01*b01*b11-a10*b01*b02),
]

function td(f)
    return total_degree(numerator(f // one(R))) + total_degree(denominator(f // one(R)))
end

println("Number of original generators: $(length(gens)). Total degrees $(map(td, gens))")
sg = simplified_generating_set(RationalFunctionField(gens))
println("Number of new generators: $(length(sg)). Total degrees $(map(td, sg))")
println(sg)
