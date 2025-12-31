using StructuralIdentifiability, RationalFunctionFields, Nemo
using Combinatorics, Random, ParamPunPam

ring, (T, a, d, dr, e, g, r, rR) = polynomial_ring(QQ, [:T, :a, :d, :dr, :e, :g, :r, :rR])

fracs = [T, T*a + T*d + T*dr + e + g - r - rR, (d*rR - dr*r)//(d - dr), (a*r - a*rR + d*e + d*g - dr*e - dr*g)//(d - dr), (a^2 + 2*a*d + d^2 + dr^2)//(a*dr + d*dr), (a*dr*r - a*dr*rR + d^2*g + d*dr*e - d*dr*g - dr^2*e)//(a*d - a*dr + d^2 - dr^2)]

if !isdefined(Main, :rff)
rff = RationalFunctionField(fracs)
end
vars = gens(rff.mqs.parent_ring_param)

total = 10000

for (i, _) in enumerate(1:total)
    i >= total && break
    V = Random.shuffle(vars)
    n = length(V)
    k = rand(1:n)
    A = vars[1:k]
    B = vars[k+1:end]
    RationalFunctionFields.groebner_basis_coeffs(rff, ordering=DegRevLex(A)*DegRevLex(B), up_to_degree=(2,2))
end

CFS = empty(fracs)
for gb in values(rff.mqs.cached_groebner_bases)
    cfs = reduce(vcat, map(f -> collect(coefficients(f)), gb))
    append!(CFS, cfs)
end

CFS = unique(CFS)
POLYS =  RationalFunctionFields.polynomial_generators(RationalFunctionField(CFS), 4)

ALL = vcat(CFS, POLYS)

simple = RationalFunctionFields.beautiful_generators(RationalFunctionField(ALL))
