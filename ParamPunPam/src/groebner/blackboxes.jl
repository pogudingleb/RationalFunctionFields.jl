
"""
    AbstractBlackboxIdeal

Blackbox ideals in `K(x)[y]` that can be evaluated at `x`.

## Interface

Subtypes of `AbstractBlackboxIdeal` must implement the following functions:

- `length(<:AbstractBlackboxIdeal)`: the number of ideal generators.
- `base_ring(<:AbstractBlackboxIdeal)`: original ground field.
- `parent(<:AbstractBlackboxIdeal)`: original polynomial ring.
- `parent_params(<:AbstractBlackboxIdeal)`: original polynomial ring of coefficients.
- `reduce_mod_p!(<:AbstractBlackboxIdeal, finite_field)`: reduces generators
  modulo `p`, where `p` is the characteristic of the given `finite_field`.
- `specialize_mod_p(<:AbstractBlackboxIdeal, point)`: specializes the ideal at
  `point` and returns the generators of the specialized ideal in Zp[y].
"""
abstract type AbstractBlackboxIdeal end

mutable struct BasicBlackboxIdeal{PolyQQX,PolyFpSmall,PolyFpLarge} <: AbstractBlackboxIdeal
    polys::Vector{PolyQQX}
    small_prime::BigInt
    polys_mod_small_p::Vector{PolyFpSmall}
    large_prime::BigInt
    polys_mod_large_p::Vector{PolyFpLarge}

    function BasicBlackboxIdeal(polys::Vector{T}) where {T}
        @assert !isempty(polys)
        Rx = parent(polys[1])
        Ra = base_ring(Rx)
        polys = filter(!iszero, polys)
        @debug "Constructing a blackbox from $(length(polys)) input polynomials"
        if Ra isa Nemo.FracField
            Ra = base_ring(Ra)
            Rlifted, _ = polynomial_ring(Ra, map(string, gens(Rx)), internal_ordering=Nemo.internal_ordering(Rx))
            polys = liftcoeffs(polys, Rlifted)
        end
        K = base_ring(Ra)
        @assert K isa Nemo.FracField{Nemo.ZZRingElem} "Only Nemo.QQ as a ground field is supported"
        new{eltype(polys), AbstractAlgebra.Generic.MPoly{Nemo.fpMPolyRingElem}, AbstractAlgebra.Generic.MPoly{Nemo.FqMPolyRingElem}}(
            polys,
            BigInt(0),
            Vector{AbstractAlgebra.Generic.MPoly{Nemo.fpMPolyRingElem}}(),
            BigInt(0),
            Vector{AbstractAlgebra.Generic.MPoly{Nemo.FqMPolyRingElem}}()
        )
    end
end

AbstractAlgebra.base_ring(ideal::BasicBlackboxIdeal) = base_ring(base_ring(first(ideal.polys)))
AbstractAlgebra.parent(ideal::BasicBlackboxIdeal) = parent(first(ideal.polys))
parent_params(ideal::BasicBlackboxIdeal) = base_ring(parent(first(ideal.polys)))
Base.length(ideal::BasicBlackboxIdeal) = length(ideal.polys)

function reduce_mod_p!(ideal::BasicBlackboxIdeal, finite_field)
    @debug "Reducing modulo $(finite_field).."
    prime = BigInt(characteristic(finite_field))
    if finite_field isa Nemo.fpField
        ideal.polys_mod_small_p = map(poly -> map_coefficients(f -> map_coefficients(c -> finite_field(c), f), poly), ideal.polys)
        ideal.small_prime = prime
    else
        @assert finite_field isa Nemo.FqField
        ideal.polys_mod_large_p = map(poly -> map_coefficients(f -> map_coefficients(c -> finite_field(c), f), poly), ideal.polys)
        ideal.large_prime = prime
    end
    nothing
end

function _specialize_mod_p(polys_mod_p, point)
    for poly in polys_mod_p
        if iszero(evaluate(leading_coefficient(poly), point))
            __throw_unlucky_cancellation()
        end
    end
    map(f -> map_coefficients(c -> evaluate(c, point), f), polys_mod_p)
end

function specialize_mod_p(ideal::BasicBlackboxIdeal, point::Vector{<:Nemo.fpFieldElem})
    finite_field = parent(first(point))
    @assert ideal.small_prime == characteristic(finite_field)
    _specialize_mod_p(ideal.polys_mod_small_p, point)
end

function specialize_mod_p(ideal::BasicBlackboxIdeal, point::Vector{<:Nemo.FqFieldElem})
    finite_field = parent(first(point))
    @assert ideal.large_prime == characteristic(finite_field)
    _specialize_mod_p(ideal.polys_mod_large_p, point)
end

function liftcoeffs(polys, newring)
    cfs = map(collect ∘ coefficients, polys)
    cfs = map(c -> map(numerator, c .* lcm(map(denominator, c))), cfs)
    newpolys = Vector{elem_type(newring)}(undef, length(polys))
    for i in 1:length(newpolys)
        G = lcm(map(denominator, collect(coefficients(polys[i]))))
        newpolys[i] = map_coefficients(c -> numerator(c * G), polys[i])
    end
    newpolys
end
