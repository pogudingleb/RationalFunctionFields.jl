using StructuralIdentifiability, RationalFunctionFields, Nemo
using Combinatorics, Random, ParamPunPam, ProgressMeter

ring, (k01, k12, k13, k14, k21, k31, k41) = polynomial_ring(QQ, [:k01, :k12, :k13, :k14, :k21, :k31, :k41])

fracs = [k01, k21 + k31 + k41, k12 + k13 + k14, k21*k31 + k21*k41 + k31*k41, k12*k13 + k12*k14 + k13*k14, k12*k31 + k12*k41 + k13*k21 + k13*k41 + k14*k21 + k14*k31, k21*k31*k41, k12*k13*k14] .// 1

if !isdefined(Main, :rff)
rff = RationalFunctionField(fracs)
end
vars = gens(rff.mqs.parent_ring_param)

total = 1000

@showprogress for (i, _) in enumerate(1:total)
    i >= total && break
    V = Random.shuffle(vars)
    n = length(V)
    k = rand(1:n)
    A = vars[1:k]
    B = vars[k+1:end]
    RationalFunctionFields.groebner_basis_coeffs(rff, ordering=DegRevLex(A)*Lex(B), up_to_degree=(5,5))
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

fields_equal(RationalFunctionField(ALL), RationalFunctionField(fracs), 0.99)
