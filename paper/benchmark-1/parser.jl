# Parsing polynomials from strings

import Nemo

function parse_polynomial_from_terms(
    ring,
    terms_exploded_str,
    str_to_var
)
    base_field = Nemo.base_ring(ring)
    n = Nemo.nvars(ring)
    coeffs = Vector{Nemo.elem_type(base_field)}(undef, length(terms_exploded_str))
    monoms = Vector{Nemo.elem_type(ring)}(undef, length(terms_exploded_str))
    for (i, t) in enumerate(terms_exploded_str)
        mon_ = one(ring)
        cf_ = one(base_field)
        for m in t
            varexp = split(m, "^")
            if length(varexp) == 1
                x = varexp[1]
                if haskey(str_to_var, x)
                    mon_ *= str_to_var[x]
                else
                    constant = parse(Rational{BigInt}, x)
                    cf_ *= base_field(constant)
                end
            else
                @assert length(varexp) == 2
                x, y = varexp
                mon_ *= str_to_var[x]^parse(Int, y)
            end
        end
        coeffs[i] = cf_
        monoms[i] = mon_
    end
    coeffs, monoms
end

function parse_poly(ring, poly_str)
    terms_str_plus = split(poly_str, '+')
    terms_str = empty(terms_str_plus)
    @assert length(terms_str_plus) > 0
    for term_str in terms_str_plus
        term_str_minus = split(term_str, '-')
        if !isempty(strip(term_str_minus[1]))
            push!(terms_str, strip(term_str_minus[1]))
        end
        append!(terms_str, map(s -> "-1*$(strip(s))", term_str_minus[2:end]))
    end
    terms_exploded_str = map(t -> map(strip, split(t, "*")), terms_str)
    str_to_var = Dict(string.(symbols(ring)) .=> gens(ring))
    cfs, exps = parse_polynomial_from_terms(
        ring,
        terms_exploded_str,
        str_to_var
    )
    m = map(f -> exponent_vector(f, 1), exps)
    ring(cfs, m)
end

begin
R, (a,b,x1,x2) = Nemo.polynomial_ring(Nemo.QQ, ["a", "b", "x1", "x2"])
@assert parse_poly(R, "0") == R(0)
@assert parse_poly(R, "242342342342423142312//3432483843848483848348384834838477") == R(242342342342423142312//3432483843848483848348384834838477)
@assert parse_poly(R, "  a  -   b     * x1  ") == a - b*x1
@assert parse_poly(R, "a*b * x1 +x2^3 - 11*x1*x2") == (a*b * x1 +x2^3 - 11*x1*x2)
end
