# Discrete logarithms in F_q^n

# Stores the field K and some information about it
# to speed up discrete logarithms
mutable struct PrecomputedField{Field, T <: Integer}
    K::Field
    ordmult::T
    factors::Vector{Pair{T, Int}}
    extensiondeg::Int

    function PrecomputedField(K::Field) where {Field}
        ordmult_big = BigInt(Nemo.order(K) - 1)
        ordmult_big > BigInt(typemax(UInt64)) &&
            throw(ArgumentError("Discrete logarithms currently support multiplicative orders up to UInt64."))
        ordmult = UInt64(ordmult_big)
        factors = map(f -> UInt64(first(f)) => Int(last(f)), collect(Primes.factor(Dict, ordmult)))
        new{Field, typeof(ordmult)}(K, ordmult, factors, Nemo.degree(K))
    end
end

# Stores preallocated buffers to speed discrete logarithms 
mutable struct DiscreteLogBuffers{Field, T <: Integer, I}
    PF::PrecomputedField{Field, T}
    xibuf::Vector{T}
    pibuf::Vector{T}
    baby::Dict{I, Int}

    function DiscreteLogBuffers(PF::PrecomputedField{Field, T}) where {Field, T}
        nfactors = length(PF.factors)
        largest = maximum(f -> first(f), PF.factors)
        baby = sizehint!(Dict{elem_type(PF.K), Int}(), Int(isqrt(BigInt(largest))))
        new{Field, T, elem_type(PF.K)}(PF, Vector{T}(undef, nfactors), Vector{T}(undef, nfactors), baby)
    end
end

# Dispatch between direct and baby-giant discrete log algorithms 
# is based on the size of the base order.
# If the order is < _threshold_direct_case(), direct algorithm is used.
_threshold_direct_case() = 2^5

# Solves a^x = y (mod p) for x
# ord is the order of a in Z/Zp, factors is an array of factors of p-1
discrete_log(a::I, y::I, buf; ord_isprime=false) where {I} =
    discrete_log(a, y, buf.PF.ordmult, buf; ord_isprime=ord_isprime)

function discrete_log(a::I, y::I, ord::T, buf; ord_isprime=false) where {I, T <: Integer}
    if ord < _threshold_direct_case()
        direct_discrete_log(a, y, ord, buf)
    else
        if !ord_isprime
            # make sure pohlig_hellman_discrete_log 
            # can not call itself recursively
            pohlig_hellman_discrete_log(a, y, ord, buf)
        else
            babystep_giantstep_discrete_log(a, y, ord, buf)
        end
    end
end

# Solves a^x = y (mod p) for x using the Pohlig-Hellman algorithm.
# 
# ord is the order of a in Z/Zp.
# factors is a dictionary of prime factors of ord (with multiplicities).
# xibuf and pibuf are buffers used to store intermediate results.
function pohlig_hellman_discrete_log(a::I, y::I, ord::T, buf) where {I, T <: Integer}
    PF = buf.PF
    @inbounds for i in 1:length(PF.factors)
        (pi, di) = PF.factors[i]
        ai, yi = a, y
        xi = zero(T)
        cc = one(ai)
        for j in 0:(di - 1)
            pij = pi^j
            cij = div(ord, pij * pi)
            aij = ai^(cij * pij)
            yij = (yi * inv(cc))^cij
            xij = discrete_log(aij, yij, pi, buf, ord_isprime=true)
            tij = xij * pij
            cc *= ai^tij
            xi += tij
        end
        buf.xibuf[i] = xi
        buf.pibuf[i] = pi^di
    end
    mod(T(Nemo.crt(map(BigInt, buf.xibuf), map(BigInt, buf.pibuf))), ord)
end

# Solves a^x = y (mod p) for x using the baby-step giant-step algorithm.
# ord is the order of a in Z/Zp.
function babystep_giantstep_discrete_log(a::I, y::I, ord::T, buf) where {I, T <: Integer}
    # the size of a giant step
    m = T(isqrt(BigInt(ord)) + 1)
    m_int = Int(m)
    baby = buf.baby
    # this does nothing if baby has enough capacity
    sizehint!(baby, m_int)
    ai = one(a)
    for i in 0:(m_int - 1)
        baby[ai] = i
        ai *= a
    end
    # find a match
    ainvm = inv(a)^m
    giantstep = y
    i = zero(T)
    while !haskey(baby, giantstep)
        giantstep *= ainvm
        i += one(T)
    end
    laststep = baby[giantstep]
    empty!(baby)
    i * m + T(laststep)
end

# Solves a^x = y (mod p).
function direct_discrete_log(a::I, y::I, ord::T, buf) where {I, T <: Integer}
    i = zero(T)
    ai = one(a)
    while ai != y
        ai *= a
        i += one(T)
    end
    i
end
