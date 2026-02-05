using Nemo
using RationalFunctionFields

R, (z1, z2, z3, z4, z5, z6, z7, z8, z9) = polynomial_ring(QQ, ["z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9"])

gens = [z5//z6, 1//z6*z4, 1//z9*z7, 1//z9*z8, (z4*z9-z6*z7)//(z5*z9-z6*z8), z6*(z1*z5*z9-z1*z6*z8-z2*z4*z9+z2*z6*z7+z3*z4*z8-z3*z5*z7)//(z5*z9-z6*z8), (z1*z5*z9-z1*z6*z8-z2*z4*z9+z2*z6*z7+z3*z4*z8-z3*z5*z7)^2//z9//(z5*z9-z6*z8)^2, z6*z9*(z4*z8-z5*z7)//(z1*z5*z9-z1*z6*z8-z2*z4*z9+z2*z6*z7+z3*z4*z8-z3*z5*z7), z6*z9*(z5*z9-z6*z8)//(z1*z5*z9-z1*z6*z8-z2*z4*z9+z2*z6*z7+z3*z4*z8-z3*z5*z7), -2*(z4*z8-z5*z7)//(z5*z9-z6*z8), -(z4*z8-z5*z7)//(z5*z9-z6*z8), 2*(z4*z9-z6*z7)//(z5*z9-z6*z8), -(z4*z8-z5*z7)^2//(z5*z9-z6*z8)^2, -(z4*z9-z6*z7)^2//(z5*z9-z6*z8)^2, 2*(z4*z9-z6*z7)*(z4*z8-z5*z7)//(z5*z9-z6*z8)^2, -z6*z9*(z4*z9-z6*z7)//(z1*z5*z9-z1*z6*z8-z2*z4*z9+z2*z6*z7+z3*z4*z8-z3*z5*z7)]

sg = simplified_generating_set(RationalFunctionField(gens))

println("Simplified generators:")
for f in sg
    println(f)
end
