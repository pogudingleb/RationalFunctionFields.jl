# ------------------------------------------------------------------------------

"""
    minimize_modp(fracs)

Takes as input a list of rational functions `fracs`.
Returns a sublist generating the same field and minimal with respect to 
inclusion. The inclusion checks are carried out modulo the provided prime `p`
(so, no explicit correctness probability guarantee)
"""
function minimize_modp(fracs, prime = 2^31 - 1)
    new_fracs = copy(fracs)
    for i in length(fracs):-1:1
        leave_one_out = vcat(new_fracs[1:(i - 1)], new_fracs[(i + 1):end])
        if isempty(leave_one_out) && !isconstant(RationalFunctionField(new_fracs))
            return new_fracs
        end
        if first(field_contains_modp(RationalFunctionField(leave_one_out, [new_fracs[i]], prime)))
            new_fracs = leave_one_out
        end
    end
    return new_fracs
end

# ------------------------------------------------------------------------------

"""
    beautiful_generators(rff::RationalFunctionField)

Given a field of rational functions `rff` returns a set of "simpler" and
standardized generators for `rff`.

Applies the following passes:
1. Filter constants,
2. Remove redundant generators.
"""
@timeit _to function beautiful_generators(
    rff::RationalFunctionField;
    discard_redundant = true,
    reversed_order = false,
    priority_variables = [],
)
    time_start = time_ns()
    fracs = dennums_to_fractions(rff.dennums)
    # Filter pass
    fracs = filter(!is_rational_func_const, fracs)
    fracs = unique(fracs)
    if isempty(fracs)
        @debug "The set of generators is empty"
        return fracs
    end
    # Remove redundant pass
    if discard_redundant
        fracs_priority = filter(f -> issubset(vars(f), priority_variables), fracs)
        fracs_rest = filter(f -> !(f in fracs_priority), fracs)
        sort!(fracs_priority, lt = rational_function_cmp)
        sort!(fracs_rest, lt = rational_function_cmp)
        fracs = vcat(fracs_priority, fracs_rest)
        @debug "The pool of fractions of size $(length(fracs))\n$(join(map(repr, fracs), ",\n"))"
        if reversed_order
            non_redundant = collect(1:length(fracs))
            for i = length(fracs):-1:1
                func = fracs[i]
                if length(non_redundant) == 1
                    continue
                end
                result = field_contains_mod_p(
                    RationalFunctionField(fracs[setdiff(non_redundant, i)]),
                    [func],
                )
                @debug "Simplification: inclusion check $func $result"
                if result[1]
                    @debug "The function $func is discarded"
                    setdiff!(non_redundant, i)
                end
            end
        else
            non_redundant = Vector{Int}()
            push!(non_redundant, 1)
            for i = 2:length(fracs)
                func = fracs[i]
                result = field_contains_mod_p(
                    RationalFunctionField(fracs[non_redundant]),
                    [func],
                )
                @debug "Simplification: inclusion check $func $result"
                if !result[1]
                    @debug "The function $func is included in the set of generators"
                    push!(non_redundant, i)
                end
            end
        end
        @debug "Out of $(length(fracs)) simplified generators there are $(length(non_redundant)) non redundant"
        fracs = fracs[non_redundant]
    end
    sort!(fracs, lt = (f, g) -> rational_function_cmp(f, g))
    spring_cleaning_pass!(fracs)
    # _runtime_logger[:id_beautifulization] += (time_ns() - time_start) / 1e9
    beautification_time = (time_ns() - time_start) / 1e9
    @debug "Generators beautified in $(beautification_time) seconds"
    return fracs
end

function spring_cleaning_pass!(fracs)
    @assert all(is_rational_func_normalized, fracs)
    for i = 1:length(fracs)
        func = fracs[i]
        num, den = unpack_fraction(func)
        if is_constant(num)
            func = den // num
        end
        num, den = unpack_fraction(func)
        if leading_coefficient(num) < 0
            func = func * leading_coefficient(num)
        end
        num, den = unpack_fraction(func)
        if !isone(leading_coefficient(num))
            func = divexact(func, leading_coefficient(num))
        end
        num, den = unpack_fraction(func)
        if is_constant(den) && is_constant(Nemo.term(num, length(num)))
            func = (num - trailing_coefficient(num)) // one(num)
        end
        fracs[i] = func
    end
    fracs
end

# ------------------------------------------------------------------------------

"""
    groebner_basis_coeffs(rff; options...)

## Options

- `ordering`: GB ordering; must be one of the orderings exported by
  `ParamPunPam` or `Groebner`.
- `up_to_degree`: a tuple of integers, the degrees of numerator and denominator.
    The result is correct up to the requested degrees.
"""
@timeit _to function groebner_basis_coeffs(
    rff::RationalFunctionField;
    seed = 42,
    ordering = Groebner.InputOrdering(),
    up_to_degree = (typemax(Int), typemax(Int)),
    rational_interpolator = :VanDerHoevenLecerf,
)
    mqs = rff.mqs
    if are_generators_zero(mqs)
        return rff
    end
    gb, fracs, new_rff = nothing, nothing, nothing
    # Check if the basis is in cache
    if haskey(mqs.cached_groebner_bases, (ordering, up_to_degree))
        @debug "Cache hit with ($ordering, $up_to_degree)!"
        gb = mqs.cached_groebner_bases[ordering, up_to_degree]
        basis_coeffs = map(collect ∘ coefficients, gb)
        fracs = collect(mapreduce(Set, union!, basis_coeffs))
        return RationalFunctionField(fracs)
    end
    #_runtime_logger[:id_calls_to_gb] += 1
    current_degrees = (2, 2)
    two_sided_inclusion = false
    while !two_sided_inclusion && all(current_degrees .<= up_to_degree)
        @debug "Computing GB with parameters up to degrees $(current_degrees)"
        runtime = @elapsed gb = ParamPunPam.paramgb(
            mqs,
            up_to_degree = current_degrees,
            ordering = ordering,
            rational_interpolator = rational_interpolator,
        )
        #_runtime_logger[:id_npoints_degree] +=
        #    ParamPunPam._runtime_data[:npoints_degree_estimation]
        #_runtime_logger[:id_npoints_interpolation] +=
        #    ParamPunPam._runtime_data[:npoints_interpolation]
        #_runtime_logger[:id_groebner_time] += runtime
        @debug "Groebner basis computed in $runtime seconds"
        basis_coeffs = map(collect ∘ coefficients, gb)
        basis_coeffs_set = mapreduce(Set, union!, basis_coeffs)
        fracs = collect(basis_coeffs_set)
        @debug "Generators up to degrees $(current_degrees) are $fracs"
        @debug "Checking two-sided inclusion modulo a prime"
        time_start = time_ns()
        new_rff = RationalFunctionField(fracs)
        # the ordering of the checks is not arbitrary - the first one can terminate earlier
        # via the algebraicity check
        if !issubfield_mod_p(rff, new_rff)
            two_sided_inclusion = false
        else
            two_sided_inclusion = issubfield_mod_p(new_rff, rff)
        end
        runtime = (time_ns() - time_start) / 1e9
        #_runtime_logger[:id_inclusion_check_mod_p] += runtime
        @debug "Inclusion checked in $(runtime) seconds. Result: $two_sided_inclusion"
        current_degrees = current_degrees .* (2, 2)
    end
    @debug "The coefficients of the Groebner basis are presented by $(length(fracs)) rational functions"
    new_rff.mqs.cached_groebner_bases[ordering, up_to_degree] = gb
    rff.mqs.cached_groebner_bases[ordering, up_to_degree] = gb
    return new_rff
end

"""
    generating_sets_fan(rff::RationalFunctionField, nbases)

Returns a set of Groebner bases for multiple different rankings of variables.

## Arguments

- `nbases`: How many bases to compute.
- Keyword `up_to_degree`: a tuple of integers, max. degrees of numerators and
  denominators. Result is correct up to the requested degrees.
"""
@timeit _to function generating_sets_fan(
    rff::RationalFunctionField{T},
    code::Integer;
    seed = 42,
    up_to_degree = (3, 3),
) where {T}
    time_start = time_ns()
    vars = gens(parent(rff.mqs))
    nbases = length(vars)
    ordering_to_generators = Dict()
    if code == 0
        return ordering_to_generators
    end
    @info "Computing $nbases Groebner bases for degrees $up_to_degree for block orderings"
    # The first basis in some ordering
    ord = InputOrdering()
    new_rff = groebner_basis_coeffs(rff, seed = seed, ordering = ord)
    cfs = beautiful_generators(new_rff)
    ordering_to_generators[ord] = cfs
    if isempty(cfs)
        return ordering_to_generators
    end
    # NOTE: maybe hide the computation of multiple bases inside
    # RationalFunctionField
    gb_rff = RationalFunctionField(cfs)
    vars = gens(parent(gb_rff.mqs))
    if length(vars) == 1
        return ordering_to_generators
    end
    if code >= 1
        for i = 1:nbases
            vars_shuffled = circshift(vars, i)
            n = length(vars_shuffled)
            # n1, n2 = div(n, 2), n - div(n, 2)
            n1, n2 = n - 1, 1
            ord = DegRevLex(vars_shuffled[1:n1]) * DegRevLex(vars_shuffled[(n1+1):end])
            @debug "Computing GB for ordering $ord"
            new_rff = groebner_basis_coeffs(
                gb_rff,
                seed = seed,
                ordering = ord,
                up_to_degree = up_to_degree,
            )
            cfs = beautiful_generators(new_rff, discard_redundant = false)
            ordering_to_generators[ord] = cfs
        end
    end
    if code >= 2
        for _ = 1:nbases
            vars_shuffled = shuffle(vars)
            n = length(vars_shuffled)
            n1, n2 = max(n - 2, 1), min(2, length(vars) - 1)
            ord = DegRevLex(vars_shuffled[1:n1]) * DegRevLex(vars_shuffled[(n1+1):end])
            @debug "Computing GB for ordering $ord"
            new_rff = groebner_basis_coeffs(
                gb_rff,
                seed = seed,
                ordering = ord,
                up_to_degree = up_to_degree,
            )
            cfs = beautiful_generators(new_rff, discard_redundant = false)
            ordering_to_generators[ord] = cfs
        end
    end
    if code >= 3
        for _ = 1:nbases
            vars_shuffled = shuffle(vars)
            n = length(vars_shuffled)
            n1 = div(n, 2)
            n2 = n - n1
            ord = DegRevLex(vars_shuffled[1:n1]) * DegRevLex(vars_shuffled[(n1+1):end])
            @debug "Computing GB for ordering $ord"
            new_rff = groebner_basis_coeffs(
                gb_rff,
                seed = seed,
                ordering = ord,
                up_to_degree = up_to_degree,
            )
            cfs = beautiful_generators(new_rff, discard_redundant = false)
            ordering_to_generators[ord] = cfs
        end
    end
    #_runtime_logger[:id_gbfan_time] = (time_ns() - time_start) / 1e9
    @info "Computed Groebner bases in $((time_ns() - time_start) / 1e9) seconds"
    return ordering_to_generators
end

SIMPLIFICATION_REGIMES = Dict(
    :weak => Dict(:poly_gens => 0, :gb_fan => 0),
    :standard => Dict(:poly_gens => 3, :gb_fan => 0),
    :strong => Dict(:poly_gens => 3, :gb_fan => 3),
)

"""
    simplified_generating_set(rff; prob_threshold = 0.99, seed = 42)

Returns a simplified set of generators for `rff`. 
Result is correct (in the Monte-Carlo sense) with probability at least `prob_threshold`.
"""
@timeit _to function simplified_generating_set(
    rff::RationalFunctionField;
    prob_threshold = 0.99,
    seed = 42,
    simplify = :standard,
    enforce_minimality = false, # being set to `true` may reduce interpretability
    check_variables = false, # almost always slows down and thus turned off
    rational_interpolator = :VanDerHoevenLecerf,
    priority_variables = [],
)
    if isconstant(rff)
        return empty([one(poly_ring(rff)) // one(poly_ring(rff))])
    end
    @info "Simplifying generating set. Simplification level: $simplify"
    #_runtime_logger[:id_groebner_time] = 0.0
    #_runtime_logger[:id_calls_to_gb] = 0
    #_runtime_logger[:id_inclusion_check_mod_p] = 0.0
    #_runtime_logger[:id_inclusion_check] = 0.0
    #_runtime_logger[:id_gbfan_time] = 0.0
    #_runtime_logger[:id_normalforms_time] = 0.0
    #_runtime_logger[:id_ranking] = 0

    # Checking membership of particular variables and adding them to the field
    if check_variables
        vars = gens(poly_ring(rff))
        containment = field_contains(rff, vars, (1.0 + prob_threshold) / 2)
        prob_threshold = (1.0 + prob_threshold) / 2
        if all(containment)
            return [v // one(poly_ring(rff)) for v in vars]
        end
        field_gens = rff.dennums
        for (v, is_contained) in zip(vars, containment)
            if is_contained
                push!(field_gens, [one(poly_ring(rff)), v])
            end
        end
        rff = RationalFunctionField(field_gens)
    end

    @assert simplify in keys(SIMPLIFICATION_REGIMES)
    simpl_params = SIMPLIFICATION_REGIMES[simplify]

    # Compute the first GB in some ordering
    new_rff = groebner_basis_coeffs(
        rff,
        seed = seed,
        rational_interpolator = rational_interpolator,
    )
    new_fracs = beautiful_generators(new_rff)
    if isempty(new_fracs)
        return new_fracs
    end

    # Compute some normal forms
    poly_generators = polynomial_generators(
        RationalFunctionField(new_fracs),
        simpl_params[:poly_gens];
        seed = seed,
    )
    append!(new_fracs, poly_generators)

    # Compute some GBs
    fan = generating_sets_fan(new_rff, simpl_params[:gb_fan]; seed = seed)
    for (ord, rff_gens) in fan
        append!(new_fracs, rff_gens)
    end

    # retaining the original generators (but not all!)
    sort!(new_fracs, lt = rational_function_cmp)
    original_fracs = generators(rff)
    sort!(original_fracs, lt = rational_function_cmp)
    merge_index = findfirst(f -> rational_function_cmp(new_fracs[end], f), original_fracs)
    if isnothing(merge_index)
        merge_index = length(original_fracs)
    else
        merge_index -= 1
    end
    @debug "Retaining $(merge_index) out of $(length(original_fracs)) original generators"
    append!(new_fracs, original_fracs[1:merge_index])
    new_fracs_unique = unique(new_fracs)
    @debug """
Final cleaning and simplification of generators. 
Out of $(length(new_fracs)) fractions $(length(new_fracs_unique)) are syntactically unique."""
    start_time = time_ns()
    runtime = @elapsed new_fracs = beautiful_generators(
        RationalFunctionField(new_fracs_unique),
        priority_variables = priority_variables,
    )
    if enforce_minimality
        before = length(new_fracs)
        new_fracs = minimize_modp(new_fracs)
        @debug "Enforcing minimality: retained $(length(new_fracs)) out of $before"
    end
    @info "Selecting generators in $((time_ns() - start_time) / 1e9)"
    @debug "Checking inclusion with probability $prob_threshold"
    runtime = @elapsed result =
        fields_equal(rff, RationalFunctionField(new_fracs), prob_threshold)
    #_runtime_logger[:id_inclusion_check] = runtime
    if !result
        @warn "Field membership check failed. Error will follow."
        throw("The new subfield generators are not correct.")
    end
    @info "Inclusion checked with probability $prob_threshold in $(runtime) seconds"
    @debug "Out of $(length(rff.mqs.nums_qq)) initial generators there are $(length(new_fracs)) independent"
    ranking = generating_set_rank(new_fracs)
    #_runtime_logger[:id_ranking] = ranking
    @debug "The ranking of the new set of generators is $ranking"
    return new_fracs
end

# ------------------------------------------------------------------------------
