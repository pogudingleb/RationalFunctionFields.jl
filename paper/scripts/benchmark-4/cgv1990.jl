using StructuralIdentifiability, RationalFunctionFields, Nemo
using Combinatorics

ring, (R, V36, V3, S, k5) = polynomial_ring(QQ, [:R, :V36, :V3, :S, :k5])

fracs = [
     R*V36 + 1//25*S*V3,
     R*V3 + S*V36,
     R*V36*k5 + 1//5*S*V36*k5,
     (V3*k5 + 5*V36*k5)//V3
] .// ring(1)

rff = RationalFunctionField(fracs)

for p in permutations(gens(rff.mqs.parent_ring_param))
    RationalFunctionFields.groebner_basis_coeffs(rff, ordering=Lex(p))
end

CFS = empty(fracs)
for gb in values(rff.mqs.cached_groebner_bases)
    cfs = reduce(vcat, map(f -> collect(coefficients(f)), gb))
    append!(CFS, cfs)
end

CFS = unique(CFS)
POLYS = filter(c -> total_degree(denominator(c)) == 0, CFS)

simple = RationalFunctionFields.beautiful_generators(RationalFunctionField(CFS))

