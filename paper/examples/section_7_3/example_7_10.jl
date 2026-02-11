using Nemo
using RationalFunctionFields
R, (z1, z2, z3, z4, z5, z6, z7) = polynomial_ring(QQ, ["z1", "z2", "z3", "z4", "z5", "z6", "z7"])

gens = [(z1*z4-z2*z3)^2//(z3^2*z7-2*z3*z4*z6+z4^2*z5), (z3^2*z7-2*z3*z4*z6+z4^2*z5)//(z1*z4-z2*z3), (z1*z3*z7-z1*z4*z6-z2*z3*z6+z2*z4*z5)//(z1*z4-z2*z3), -(z1^2*z7-2*z1*z2*z6+z2^2*z5)//(z1*z4-z2*z3), -(z1^2*z7-2*z1*z2*z6+z2^2*z5)//(z3^2*z7-2*z3*z4*z6+z4^2*z5), -(z1*z3*z7-z1*z4*z6-z2*z3*z6+z2*z4*z5)//(z1*z4-z2*z3), 2*(z1*z3*z7-z1*z4*z6-z2*z3*z6+z2*z4*z5)//(z3^2*z7-2*z3*z4*z6+z4^2*z5), -(z1*z4-z2*z3)*(z1*z3*z7-z1*z4*z6-z2*z3*z6+z2*z4*z5)//(z3^2*z7-2*z3*z4*z6+z4^2*z5), -z5*z7+z6^2, -z1*z4+z2*z3]

rff = RationalFunctionField(gens)
sg = simplified_generating_set(rff, enforce_minimality = true)

println("Simplified generators:")
for f in sg
    println(f)
end
