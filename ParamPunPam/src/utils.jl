
# Keep track of some statistics
const _runtime_data = Dict()

@noinline __throw_unlucky_cancellation() =
    throw(AssertionError("Unlucky cancellation of coefficients in the Groebner basis!"))

@noinline __throw_something_went_wrong(msg) =
    throw(AssertionError("Something went wrong when computing Groebner bases.\n$msg"))

# Progress bars
const _progressbar_color = :cyan
const _progressbar_value_color = :cyan # :light_grey
const _progressbar_spinner = "⌜⌝⌟⌞"
const _is_progressbar_enabled_globally = Ref{Bool}(true)
enable_progressbar(flag::Bool) = _is_progressbar_enabled_globally[] = flag
is_progressbar_enabled() =
    _is_progressbar_enabled_globally[] && Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

# Some other utils..
function evaluate_frac(f, x)
    n, d = numerator(f), denominator(f)
    isone(d) && return evaluate(n, x)
    evaluate(n, x) // evaluate(d, x)
end

lift_modular_elem(c) = BigInt(lift(Nemo.ZZ, c))

function crt_bigint(remainders::Vector{BigInt}, moduli::Vector{BigInt})
    @assert !isempty(remainders)
    @assert length(remainders) == length(moduli)
    x = mod(remainders[1], moduli[1])
    modulus = moduli[1]
    for i in 2:length(remainders)
        mi = moduli[i]
        delta = mod(remainders[i] - mod(x, mi), mi)
        step = mod(delta * invmod(mod(modulus, mi), mi), mi)
        x += modulus * step
        modulus *= mi
    end
    mod(x, modulus)
end
