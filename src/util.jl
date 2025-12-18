
"""
    dennums_to_fractions(dennums)
    
Returns the field generators represented by fractions.

Input: an array of arrays of polynomials, as in 
`[[f1, f2, f3, ...], [g1, g2, g3, ...], ...]`

Output: an array of fractions
`[f2/f1, f3/f1, ..., g2/g1, g3/g1, ...]`
"""
function dennums_to_fractions(dennums::Vector{Vector{T}}) where {T}
    fractions = Vector{AbstractAlgebra.Generic.FracFieldElem{T}}()
    for dni in dennums
        #pivot_ind = findmin(p -> (total_degree(p), length(p)), filter(!iszero, dni))[2]
        #pivot = dni[pivot_ind]
        append!(fractions, [c // dni[1] for c in dni[2:end]])
    end
    return fractions
end

# ------------------------------------------------------------------------------

"""
    fractions_to_dennums(fractions)
    
Returns the field generators represented by lists of denominators and
numerators.

Input: an array of fractions, as in
`[f2/f1, f3/f1, ..., g2/g1, g3/g1, ...]`

Output: an array of arrays of polynomials,
`[[f1, f2, f3, ...], [g1, g2, g3, ...], ...]`
"""
function fractions_to_dennums(fractions)
    return map(f -> [denominator(f), numerator(f)], fractions)
end

# ------------------------------------------------------------------------------

"""
    merge_results(outer::Vector{Bool}, inner::Vector{Bool})

Returns a list `res` of Bools of the length as outer such that `res[i]` is true
iff `outer[i]` is true and `inner[j]` is true where `j` is the index of `outer[i]`
among the true values in `outer`
"""

function merge_results(outer::Vector{Bool}, inner::Vector{Bool})
    @assert length(inner) == count(outer)
    result = copy(outer)
    inner_index = 1

    for i = 1:length(result)
        if result[i]
            if !inner[inner_index]
                result[i] = false
            end
            inner_index += 1
        end
    end
    return result
end

# ------------------------------------------------------------------------------

# Feels like inventing a bicycle
function squarefree_part(p)
    g = reduce(gcd, map(x -> derivative(p, x), vars(p)), init = p)
    return divexact(p, g)
end

# ------------------------------------------------------------------------------

function cancel_gcds(polys::Vector)
    cancelled_polys = [squarefree_part(p) for p in polys]
    for (i, p) in enumerate(cancelled_polys)
        for j = (i+1):length(cancelled_polys)
            cancelled_polys[j] = divexact(cancelled_polys[j], gcd(cancelled_polys[j], p))
        end
    end
    @debug "Degrees before taking product $(map(total_degree, cancelled_polys))"
    @debug "Length before taking product $(map(length, cancelled_polys))"
    return filter(p -> total_degree(p) > 0, cancelled_polys)
end

# ------------------------------------------------------------------------------

function insert_at_indices(arr::Vector, indices::Vector{Int}, elem)
    result = empty(arr)
    idx_arr = 1
    for i = 1:(length(arr)+length(indices))
        if i in indices
            push!(result, elem)
        else
            push!(result, arr[idx_arr])
            idx_arr += 1
        end
    end
    return result
end

# ------------------------------------------------------------------------------

# a Groebner basis of an ideal in the parameter ring.
# Assumes that there is no division by zero
function normalize_coefficients(poly, coeff_relations)
    res = zero(parent(poly))
    for (c, m) in zip(coefficients(poly), monomials(poly))
        c_num = numerator(c)
        _, c_num = divrem(c_num, coeff_relations)
        res += c_num // denominator(c) * m
    end
    return res
end

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

# ------------------------------------------------------------------------------

"""
    parent_ring_change(poly, new_ring)

Converts a polynomial to a different polynomial ring
Input
- `poly` - a polynomial to be converted
- `new_ring` - a polynomial ring such that every variable name appearing in poly appears among the generators

Output:
- a polynomial in `new_ring` “equal” to `poly`
"""
function parent_ring_change(
    poly::MPolyRingElem,
    new_ring::MPolyRing;
    matching = :byname,
    shift = 0,
)
    old_ring = parent(poly)
    # Construct a mapping for the variable indices.
    # Zero indicates no image of the old variable in the new ring  
    var_mapping = zeros(Int, max(nvars(old_ring), nvars(new_ring)))
    if matching === :byname
        old_symbols, new_symbols = symbols(old_ring), symbols(new_ring)
        for i = 1:length(old_symbols)
            u = old_symbols[i]
            found = findfirst(v -> (u === v), new_symbols)
            isnothing(found) && continue
            var_mapping[i] = found
        end
    elseif matching === :byindex
        var_mapping[1:(nvars(new_ring)-shift)] .= (1+shift):nvars(new_ring)
    else
        throw(Base.ArgumentError("Unknown matching type: $matching"))
    end
    # Hoist the compatibility check out of the loop
    for i = 1:nvars(old_ring)
        if degree(poly, i) > 0 && iszero(var_mapping[i])
            throw(
                Base.ArgumentError(
                    """
                    The polynomial $poly contains a variable $(gens(old_ring)[i]) not present in the new ring.
                    New ring variables are $(gens(new_ring)))""",
                ),
            )
        end
    end
    bring = base_ring(new_ring)
    exps = Vector{Vector{Int}}(undef, length(poly))
    coefs = map(c -> bring(c), coefficients(poly))
    @inbounds for i = 1:length(poly)
        evec = exponent_vector(poly, i)
        new_exp = zeros(Int, nvars(new_ring))
        for i = 1:length(evec)
            iszero(var_mapping[i]) && continue
            new_exp[var_mapping[i]] = evec[i]
        end
        exps[i] = new_exp
    end
    return new_ring(coefs, exps)
end

function parent_ring_change(
    f::Generic.FracFieldElem{<:MPolyRingElem},
    new_ring::MPolyRing;
    matching = :byname,
)
    n, d = unpack_fraction(f)
    return parent_ring_change(n, new_ring; matching = matching) //
           parent_ring_change(d, new_ring; matching = matching)
end

# ------------------------------------------------------------------------------

function total_degree_frac(f::Generic.FracFieldElem{<:MPolyRingElem})
    return sum(map(total_degree, unpack_fraction(f)))
end

function total_degree_frac(f::MPolyRingElem)
    return total_degree(f)
end

# ------------------------------------------------------------------------------

"""
    compare_rational_func_by(f, g, by)

Returns 
- `-1` if `f < g`,
- `0` if `f = g`, and 
- `1` if `f > g`.

Functions' numerators and denominators are compared using `by`.
"""
function compare_rational_func_by(f, g, by::By, priority = :numerator) where {By}
    # Specializes on the type of `by`
    numf, denf = unpack_fraction(f)
    numg, deng = unpack_fraction(g)
    keynumf, keydenf = by(numf), by(denf)
    keynumg, keydeng = by(numg), by(deng)
    if priority === :numerator
        keynumf < keynumg && return -1
        keynumf > keynumg && return 1
        keydenf < keydeng && return -1
        keydenf > keydeng && return 1
    elseif priority === :denominator
        keydenf < keydeng && return -1
        keydenf > keydeng && return 1
        keynumf < keynumg && return -1
        keynumf > keynumg && return 1
    elseif priority === :additive
        keydenf + keynumf < keynumg + keydeng && return -1
        keydenf + keynumf > keynumg + keydeng && return 1
    else
        throw(DomainError("Unknown value for keyword argument priority", priority))
    end
    return 0
end

# ------------------------------------------------------------------------------

"""
    select_pivots(M::MatElem)

Takes as input a matrix M in the reduced row echelon form and returns
tow lists: of the pivot indices and the non-pivot ones
"""
function select_pivots(M::MatElem)
    @assert is_rref(M)
    j = 1
    (nrows, ncols) = size(M)
    nonpivots = Vector{Int}()
    pivots = Vector{Int}()
    for i = 1:ncols
        pivot = false
        for k = j:nrows
            if !iszero(M[k, i])
                j = k + 1
                pivot = true
            end
        end
        if pivot
            push!(pivots, i)
        else
            push!(nonpivots, i)
        end
    end
    return pivots, nonpivots
end

# -----------------------------------------------------------------------------


"""
    jacobian(ratfuncs, point)

Computes the evaluation of the jacobian of `ratfuncs` at point `point`
"""
function jacobian(ratfuncs, point)
    parent_polyring = parent(numerator(first(ratfuncs)))
    F = base_ring(parent_polyring)
    base_vars = gens(parent_polyring)
    @assert length(base_vars) == length(point)
    S = matrix_space(F, length(base_vars), length(ratfuncs))
    J = zero(S)
    for (i, f) in enumerate(ratfuncs)
        for (j, x) in enumerate(base_vars)
            J[j, i] = evaluate(derivative(f, x), point)
        end
    end
    return J
end

# ------------------------------------------------------------------------------

"""
    gen_tag_name(base; stop_words)
    gen_tag_names(n, base; stop_words)

Generates a string which will not collide with the words in `stop_words`.

## Arguments

- `n`: Generates a sequence of unique strings of length `n`
- `base`: A string or a vector of strings, the base for the generated sequence
- `stop_words`: A vector of strings, stop words
"""
function gen_tag_name(base = "T"; stop_words = Vector{String}())
    return first(gen_tag_names(1, base, stop_words = stop_words))
end

function gen_tag_names(n::Integer, base = "T"; stop_words = Vector{String}())
    sequence = if base isa Vector{String}
        @assert n == length(base)
        base
    else
        repeat([base], n)
    end
    while true
        rand_token = Int(rand(UInt8))
        sequence = map(c -> "$(rand_token)__$c", sequence)
        sequence = map(ic -> "$(ic[2])_$(ic[1])", enumerate(sequence))
        if all(elem -> !(elem in stop_words), sequence)
            break
        end
    end
    return sequence
end

# ------------------------------------------------------------------------------

function str_to_var(s::String, ring::MPolyRing)
    ind = findfirst(v -> (string(v) == s), symbols(ring))
    if ind === nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end

# ------------------------------------------------------------------------------

function var_to_str(v::MPolyRingElem; xs = gens(parent(v)))
    ind = findfirst(vv -> vv == v, xs)
    return string(symbols(parent(v))[ind])
end

# ------------------------------------------------------------------------------

"""
    _reduce_mod_p(f, p)

Reduces a polynomial/rational function over Q modulo p
"""
function _reduce_mod_p(poly::QQMPolyRingElem, p::Int)
    den = denominator(poly)
    num = change_base_ring(Nemo.ZZ, den * poly)
    if Nemo.Native.GF(p)(den) == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
    end
    return change_base_ring(Nemo.Native.GF(p), num) * (1 // Nemo.Native.GF(p)(den))
end

function _reduce_mod_p(rat::Generic.FracFieldElem{QQMPolyRingElem}, p::Int)
    num, den = map(poly -> _reduce_mod_p(poly, p), [numerator(rat), denominator(rat)])
    if den == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $rat"))
    end
    return num // den
end

#--------------------------------------

function is_rational_func_const(f)
    is_constant(numerator(f)) && is_constant(denominator(f))
end

function is_rational_func_normalized(f)
    leading_coefficient(denominator(f)) > 0 && isone(gcd(numerator(f), denominator(f)))
end
