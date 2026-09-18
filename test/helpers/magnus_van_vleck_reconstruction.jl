using FloquetExpansions

include(joinpath(@__DIR__, "magnus_log_recurrence.jl"))

const FE_MVVR = FloquetExpansions

function magnus_van_vleck_reconstruction(
  plan::FE_MVVR.BlochProjectionPlan{H},
  bloch::FE_MVVR.BlochProjectionResult{H,T};
  product,
  zero_component::T,
  simplifier=identity,
) where {H,T}
  order = length(bloch.effective) - 1
  length(bloch.wave) == order ||
    throw(ArgumentError("Bloch wave/effective truncations are inconsistent"))
  plan.order == order + 1 ||
    throw(ArgumentError("Bloch plan/result truncations are inconsistent"))

  identity_component = one(first(bloch.effective))
  static_factor = T[identity_component]
  normalized_embedding = Vector{Dict{H,T}}(undef, order)
  log_embedding = Vector{Dict{H,T}}(undef, order)
  logarithmic_derivative = Vector{Dict{H,T}}(undef, order)
  ad_coefficients = [Vector{Dict{H,T}}(undef, order) for _ in 1:max(order - 1, 0)]
  bernoulli = mlr_bernoulli_numbers(order)
  counts = FE_MVVR.BlochVanVleckCounts()

  for grade in 1:order
    n = grade - 1

    # U = (1 + Y) N.  At this point N_grade is still unknown.
    normalized_n = copy(bloch.wave[grade])
    for j in 1:(grade - 1)
      counts.factor_products += 1
      correction = FE_MVVR.bloch_vv_right_static_product(
        bloch.wave[j], static_factor[grade - j + 1], product, simplifier, counts
      )
      normalized_n = FE_MVVR.bloch_vv_add(normalized_n, correction, simplifier)
    end

    # Right logarithmic derivative A U = U'.  The current static N_grade
    # enters A_{grade-1} linearly and does not enter any nested-ad term below.
    derivative_n = FE_MVVR.bloch_vv_scale(grade, normalized_n, simplifier)
    for j in 1:n
      correction = FE_MVVR.bloch_vv_periodic_product(
        logarithmic_derivative[n - j + 1],
        normalized_embedding[j],
        product,
        simplifier,
        counts,
      )
      derivative_n = FE_MVVR.bloch_vv_add(
        derivative_n, FE_MVVR.bloch_vv_scale(-1, correction, simplifier), simplifier
      )
    end

    omega_derivative_n = copy(derivative_n)
    for depth in 1:n
      coefficient = Dict{H,T}()
      for j in 1:(n - depth + 1)
        right =
          depth == 1 ? logarithmic_derivative[n - j + 1] :
          ad_coefficients[depth - 1][n - j + 1]
        contribution = mlr_periodic_commutator(
          log_embedding[j], right, product, simplifier, counts
        )
        coefficient = FE_MVVR.bloch_vv_add(coefficient, contribution, simplifier)
      end
      ad_coefficients[depth][n + 1] = coefficient
      weight = bernoulli[depth + 1] / factorial(depth)
      omega_derivative_n = FE_MVVR.bloch_vv_add(
        omega_derivative_n,
        FE_MVVR.bloch_vv_scale(weight, coefficient, simplifier),
        simplifier,
      )
    end

    log_n = FE_MVVR.bloch_vv_scale(1 // grade, omega_derivative_n, simplifier)
    static_n = simplifier(-get(log_n, plan.zero_harmonic, zero_component))::T
    push!(static_factor, static_n)

    if !iszero(static_n)
      FE_MVVR.bloch_vv_accumulate!(normalized_n, plan.zero_harmonic, static_n)
      normalized_n = FE_MVVR.bloch_vv_simplify_embedding(normalized_n, simplifier)
      derivative_n = FE_MVVR.bloch_vv_add(
        derivative_n,
        Dict(plan.zero_harmonic => simplifier(grade * static_n)::T),
        simplifier,
      )
      FE_MVVR.bloch_vv_accumulate!(log_n, plan.zero_harmonic, static_n)
      log_n = FE_MVVR.bloch_vv_simplify_embedding(log_n, simplifier)
    end

    normalized_embedding[grade] = normalized_n
    logarithmic_derivative[grade] = derivative_n
    log_embedding[grade] = log_n
  end

  inverse_static_factor = FE_MVVR.bloch_vv_static_series_inverse(
    static_factor, order, product, counts; simplifier
  )
  right_transformed = FE_MVVR.bloch_vv_static_series_product(
    bloch.effective, static_factor, order, product, counts; simplifier
  )
  effective = FE_MVVR.bloch_vv_static_series_product(
    inverse_static_factor, right_transformed, order, product, counts; simplifier
  )

  return FE_MVVR.BlochVanVleckResult(
    static_factor,
    inverse_static_factor,
    normalized_embedding,
    log_embedding,
    effective,
    counts,
  )
end
