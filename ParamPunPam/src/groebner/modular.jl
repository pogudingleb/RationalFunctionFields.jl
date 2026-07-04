
# Rateher smooth prime (Goldilocks)
# p - 1 = 2^32 * 3 * 5 * 17 * 257 * 65537
const _default_modular_prime = big(2)^64 - big(2)^32 + 1
const _large_interpolation_prime = Primes.nextprime(big(2)^255)

finite_field_from_prime(prime::Integer) =
    prime <= typemax(UInt64) ? Nemo.Native.GF(UInt(prime)) : Nemo.GF(BigInt(prime))

requested_polynomial_interpolator(polynomial_interpolator) = begin
    polynomial_interpolator in (:PrimesBenOrTiwari, PrimesBenOrTiwari) && return PrimesBenOrTiwari
    polynomial_interpolator in (:KronBenOrTiwari, KronBenOrTiwari) && return KronBenOrTiwari
    error("Unknown polynomial interpolator: $(polynomial_interpolator)")
end

is_polynomial_interpolation_feasible(::Type{PrimesBenOrTiwari}, K, Ds::Vector{<:Integer}) =
    is_interpolation_feasible(maximum(Ds), K, length(Ds))
is_polynomial_interpolation_feasible(::Type{KronBenOrTiwari}, K, Ds::Vector{<:Integer}) =
    is_kron_interpolation_feasible(Ds, K)

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

function maybe_promote_interpolation_prime!(
    mt::ModularTracker,
    interpolation_prime_bits,
    polynomial_interpolator,
    interpolation_degrees::Vector{<:Integer}
)
    selected_interpolator = requested_polynomial_interpolator(polynomial_interpolator)
    if selected_interpolator === KronBenOrTiwari
        return false
    end
    current_prime = BigInt(characteristic(mt.finite_field))
    if interpolation_prime_bits === :auto
        current_prime > typemax(UInt64) && return false
        is_polynomial_interpolation_feasible(selected_interpolator, mt.finite_field, interpolation_degrees) && return false
        mt.finite_field = finite_field_from_prime(_large_interpolation_prime)
        return true
    end
    @assert interpolation_prime_bits in (64, 256)
    is_polynomial_interpolation_feasible(selected_interpolator, mt.finite_field, interpolation_degrees) && return false
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
