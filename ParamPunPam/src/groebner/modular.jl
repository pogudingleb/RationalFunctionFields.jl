
const _default_modular_prime = big(2)^64 - 59
const _large_interpolation_prime = Primes.nextprime(big(2)^255)

finite_field_from_prime(prime::Integer) =
    prime <= typemax(UInt64) ? Nemo.Native.GF(UInt(prime)) : Nemo.GF(BigInt(prime))

mutable struct ModularTracker
    # Current finite field
    finite_field::Any
    # The product of all used prime numbers but the last one
    modulo::BigInt

    used_primes::Vector{BigInt}

    function ModularTracker(blackbox)
        finite_field = finite_field_from_prime(_default_modular_prime)
        new(finite_field, BigInt(1), BigInt[])
    end
end

function maybe_promote_interpolation_prime!(mt::ModularTracker, interpolation_prime_bits, D::Integer, n::Integer)
    current_prime = BigInt(characteristic(mt.finite_field))
    if interpolation_prime_bits === :auto
        current_prime > typemax(UInt64) && return false
        is_interpolation_feasible(D, mt.finite_field, n) && return false
        mt.finite_field = finite_field_from_prime(_large_interpolation_prime)
        return true
    end
    @assert interpolation_prime_bits in (64, 256)
    if interpolation_prime_bits == 64
        return false
    end
    if current_prime <= typemax(UInt64)
        mt.finite_field = finite_field_from_prime(_large_interpolation_prime)
        return true
    end
    false
end

function find_next_lucky_prime!(mt::ModularTracker)
    current_prime = BigInt(characteristic(mt.finite_field))
    next_prime = Primes.prevprime(current_prime - 1)
    push!(mt.used_primes, current_prime)
    mt.finite_field = finite_field_from_prime(next_prime)
    nothing
end
