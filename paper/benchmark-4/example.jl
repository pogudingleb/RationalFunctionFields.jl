using RationalFunctionFields, ParamPunPam, Nemo

R, (a,b,c,A,B,C) = polynomial_ring(QQ, [:a,:b,:c,:A,:B,:C])

gens = [
    a^2 + b^2 - 2*a*b*A,
    b^2 + c^2 - 2*b*c*B,
    c^2 + a^2 - 2*c*a*C,
]

rff = RationalFunctionField(gens)

rff_gb = RationalFunctionFields.groebner_basis_coeffs(rff, up_to_degree=(2,2))
fracs = RationalFunctionFields.dennums_to_fractions(rff_gb.dennums)

sim = RationalFunctionFields.simplified_generating_set(rff, simplify=:strong)

R, X = polynomial_ring(QQ, :X => (1:4, 1:4))

f = Nemo.charpoly(matrix(X))
cfs = collect(coefficients(f))
