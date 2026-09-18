using FloquetExpansions

const FE_MLR = FloquetExpansions

function mlr_bernoulli_numbers(order::Int)
  order >= 0 || throw(ArgumentError("order must be nonnegative"))
  values = Rational{Int}[1 // 1]
  for n in 1:order
    total = sum(binomial(n + 1, k) * values[k + 1] for k in 0:(n - 1); init=0 // 1)
    push!(values, -total / (n + 1))
  end
  return values
end

function mlr_periodic_commutator(left, right, product, simplifier, counts)
  forward = FE_MLR.bloch_vv_periodic_product(left, right, product, simplifier, counts)
  backward = FE_MLR.bloch_vv_periodic_product(right, left, product, simplifier, counts)
  return FE_MLR.bloch_vv_add(
    forward, FE_MLR.bloch_vv_scale(-1, backward, simplifier), simplifier
  )
end

function connected_log_magnus_recurrence(
  normalized::Vector{Dict{H,T}};
  product,
  zero_component::T,
  simplifier=identity,
) where {H,T}
  order = length(normalized)
  order == 0 && return Vector{Dict{H,T}}(), FE_MLR.BlochVanVleckCounts()

  counts = FE_MLR.BlochVanVleckCounts()
  logarithmic_derivative = Vector{Dict{H,T}}(undef, order)

  # A U = U', where U = 1 + sum_{n>=1} lambda^n U_n.
  for n in 0:(order - 1)
    coefficient = FE_MLR.bloch_vv_scale(n + 1, normalized[n + 1], simplifier)
    for j in 1:n
      correction = FE_MLR.bloch_vv_periodic_product(
        logarithmic_derivative[n - j + 1], normalized[j], product, simplifier, counts
      )
      coefficient = FE_MLR.bloch_vv_add(
        coefficient, FE_MLR.bloch_vv_scale(-1, correction, simplifier), simplifier
      )
    end
    logarithmic_derivative[n + 1] = coefficient
  end

  bernoulli = mlr_bernoulli_numbers(order)
  omega = Vector{Dict{H,T}}(undef, order)
  ad_coefficients = [Vector{Dict{H,T}}(undef, order) for _ in 1:max(order - 1, 0)]

  for n in 0:(order - 1)
    derivative = copy(logarithmic_derivative[n + 1])

    for depth in 1:n
      coefficient = Dict{H,T}()
      for j in 1:(n - depth + 1)
        right = depth == 1 ? logarithmic_derivative[n - j + 1] : ad_coefficients[depth - 1][n - j + 1]
        contribution = mlr_periodic_commutator(
          omega[j], right, product, simplifier, counts
        )
        coefficient = FE_MLR.bloch_vv_add(coefficient, contribution, simplifier)
      end
      ad_coefficients[depth][n + 1] = coefficient
      weight = bernoulli[depth + 1] / factorial(depth)
      derivative = FE_MLR.bloch_vv_add(
        derivative, FE_MLR.bloch_vv_scale(weight, coefficient, simplifier), simplifier
      )
    end

    omega[n + 1] = FE_MLR.bloch_vv_scale(1 // (n + 1), derivative, simplifier)
  end

  return omega, counts
end
