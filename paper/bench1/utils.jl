# ------------------------------------------------------------------------------
# Utilities from SI.jl

"""
    eval_at_dict(f, d)

Evaluates a polynomial/rational function on a dictionary of type `var => val` and missing values are replaced with zeroes
"""
function eval_at_dict(poly::P, d::Dict{P,<:RingElem}) where {P<:MPolyRingElem}
    R = parent(first(values(d)))
    point = [get(d, v, zero(R)) for v in gens(parent(poly))]
    return evaluate(poly, point)
end

function eval_at_dict(poly::P, d::Dict{P,S}) where {P<:MPolyRingElem,S<:Real}
    R = parent(poly)
    @assert R == parent(first(keys(d)))
    xs = gens(parent(first(keys(d))))
    xs_sym = [get(d, x, 0.0) for x in xs if string(x) in map(string, gens(R))]
    accum = zero(valtype(d))
    for t in terms(poly)
        cf = coeff(t, 1)
        # Nemo.QQ --> Rational{BigInt}
        # NOTE: what about Nemo.GF|ZZ?
        cf_ = BigInt(numerator(cf)) // BigInt(denominator(cf))
        ex = exponent_vector(t, 1)
        accum += cf_ * prod(xs_sym .^ ex)
    end
    return accum
end

function eval_at_dict(
    rational::Generic.FracFieldElem{T},
    d::Dict{T,V},
) where {T<:MPolyRingElem,V}
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) / eval_at_dict(g, d)
end

function eval_at_dict(
    rational::Generic.FracFieldElem{<:T},
    d::Dict{T,<:RingElem},
) where {T<:MPolyRingElem}
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) * inv(eval_at_dict(g, d))
end

function eval_at_dict(
    rational::Generic.FracFieldElem{<:P},
    d::Dict{<:P,<:Union{<:Generic.FracFieldElem,<:P}},
) where {P<:MPolyRingElem}
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) // eval_at_dict(g, d)
end

# ------------------------------------------------------------------------------

function unpack_fraction(f::MPolyRingElem)
    return (f, one(parent(f)))
end

function unpack_fraction(f::Generic.FracFieldElem{<:MPolyRingElem})
    return (numerator(f), denominator(f))
end
