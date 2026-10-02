const Harmonics{X} = Dict{Int,X}

struct NativeRecurrence{X,O,T}
  E::Vector{X}
  S::Vector{X}
  slots::Vector{X}
  H::Vector{O}
  kick::Vector{Harmonics{X}}
  channels::Vector{GradedChannel{T}}
  known_gram::Vector{Matrix{T}}
end

function harmonic_combine(A::Harmonics{X}, B::Harmonics{X}, α::Number, β::Number) where {X}
  C = Harmonics{X}()
  for (k, value) in A
    C[k] = α * value
  end
  for (k, value) in B
    C[k] = haskey(C, k) ? C[k] + β * value : β * value
  end
  return C
end

function harmonic_product(A::Harmonics{X}, B::Harmonics{X}) where {X}
  C = Harmonics{X}()
  for (ka, a) in A, (kb, b) in B
    k = ka + kb
    C[k] = haskey(C, k) ? C[k] + a * b : a * b
  end
  return C
end

harmonic_average(A::Harmonics{X}, zero_element::X) where {X} = get(A, 0, zero_element)

function harmonic_integral(A::Harmonics{X}) where {X}
  return Harmonics{X}(k => (im // k) * value for (k, value) in A if k != 0)
end

function harmonic_derivative(A::Harmonics{X}) where {X}
  return Harmonics{X}(k => (-im * k) * value for (k, value) in A if k != 0)
end

function harmonic_right(A::Harmonics{X}, B::X) where {X}
  return Harmonics{X}(k => value * B for (k, value) in A)
end

function harmonic_with_static(A::Harmonics{X}, value::X) where {X}
  result = copy(A)
  result[0] = value
  return result
end

kick_zero(::Type{X}, N::Int) where {X} = [Harmonics{X}() for _ in 0:N]

function kick_identity(identity::X, N::Int) where {X}
  result = kick_zero(X, N)
  result[1] = Harmonics{X}(0 => identity)
  return result
end

function kick_add(
  A::Vector{Harmonics{X}}, B::Vector{Harmonics{X}}, N::Int, α::Number, β::Number
) where {X}
  return [harmonic_combine(A[k], B[k], α, β) for k in 1:(N + 1)]
end

function kick_scale(A::Vector{Harmonics{X}}, N::Int, α::Number) where {X}
  return [Harmonics{X}(k => α * value for (k, value) in A[n]) for n in 1:(N + 1)]
end

function kick_mul(A::Vector{Harmonics{X}}, B::Vector{Harmonics{X}}, N::Int) where {X}
  out = kick_zero(X, N)
  for n in 0:N
    coefficient = Harmonics{X}()
    for j in 0:n
      j + 1 <= length(A) || continue
      n - j + 1 <= length(B) || continue
      product = harmonic_product(A[j + 1], B[n - j + 1])
      coefficient = harmonic_combine(coefficient, product, 1, 1)
    end
    out[n + 1] = coefficient
  end
  return out
end

function kick_commutator(A::Vector{Harmonics{X}}, B::Vector{Harmonics{X}}, N::Int) where {X}
  return kick_add(kick_mul(A, B, N), kick_mul(B, A, N), N, 1, -1)
end

function kick_exp(A::Vector{Harmonics{X}}, identity::X, N::Int) where {X}
  result = kick_identity(identity, N)
  power = kick_identity(identity, N)
  for p in 1:N
    power = kick_mul(power, A, N)
    result = kick_add(result, kick_scale(power, N, 1 // factorial(p)), N, 1, 1)
  end
  return result
end

function kick_log(K::Vector{Harmonics{X}}, identity::X, N::Int) where {X}
  A = kick_add(K, kick_identity(identity, N), N, 1, -1)
  result = kick_zero(X, N)
  power = A
  for p in 1:N
    result = kick_add(result, kick_scale(power, N, (-1)^(p + 1) // p), N, 1, 1)
    p == N || (power = kick_mul(power, A, N))
  end
  return result
end

function static_exp(S::Vector{X}, identity::X, N::Int, sign::Int) where {X}
  A = kick_zero(X, N)
  for n in 1:min(N, length(S))
    A[n + 1] = Harmonics{X}(0 => sign * S[n])
  end
  return kick_exp(A, identity, N)
end

function intrinsic_offset(
  kick::Vector{Harmonics{X}}, S::Vector{X}, order::Int, identity::X, zero_element::X
) where {X}
  K = kick_mul(kick, static_exp(S, identity, order, -1), order)
  return harmonic_average(kick_log(K, identity, order)[order + 1], zero_element)
end

function bf_truncated_kick(Y::Vector{Harmonics{X}}, order::Int) where {X}
  truncated = [copy(Y[k]) for k in 1:(order + 1)]
  delete!(truncated[order + 1], 0)
  return truncated
end

function bf_static_residual(
  L::Harmonics{X}, Y::Vector{Harmonics{X}}, E::Vector{X}, order::Int, Yn::Harmonics{X}
) where {X}
  residual = harmonic_product(L, Yn)
  for j in 1:order
    Yj = j == order ? Yn : Y[j + 1]
    residual = harmonic_combine(residual, harmonic_right(Yj, E[order - j + 1]), 1, -1)
  end
  return residual
end

function hd_truncated_generator(G::Vector{Harmonics{X}}, order::Int) where {X}
  truncated = kick_zero(X, order)
  for n in 1:order
    truncated[n + 1] = copy(G[n + 1])
  end
  delete!(truncated[order + 1], 0)
  return truncated
end

function hd_generator_series(L::Harmonics{X}, G::Vector{Harmonics{X}}, N::Int) where {X}
  generator = kick_zero(X, N)
  for n in 1:min(N, length(G) - 1)
    generator[n + 1] = G[n + 1]
  end
  Lseries = kick_zero(X, N)
  Lseries[1] = L

  effective = [copy(value) for value in Lseries]
  term = Lseries
  for k in 1:N
    term = kick_commutator(generator, term, N)
    effective = kick_add(effective, kick_scale(term, N, (-1)^k // factorial(k)), N, 1, 1)
  end

  drift = kick_zero(X, N)
  for n in 0:N
    n + 2 <= length(G) || continue
    drift[n + 1] = harmonic_derivative(G[n + 2])
  end
  term = drift
  effective = kick_add(effective, term, N, 1, -1)
  for k in 1:N
    term = kick_commutator(generator, term, N)
    effective = kick_add(
      effective, kick_scale(term, N, -(-1)^k // factorial(k + 1)), N, 1, 1
    )
  end
  return effective
end

function initial_native_channels(
  leading::AbstractMatrix{T}, weights::AbstractVector
) where {T<:Number}
  length(weights) == size(leading, 2) ||
    throw(DimensionMismatch("one rate weight is required per leading channel"))
  return GradedChannel{T}[
    GradedChannel{T}(0, T(weights[j]), [Vector{T}(leading[:, j])]) for j in axes(leading, 2)
  ]
end

function check_native_coefficient(representation, E, channels, order, dimension, tolerance)
  reconstructed = gram_coefficient(channels, order, dimension)
  native_matches(
    representation, native_kossakowski(representation, E), reconstructed, 20 * tolerance
  ) || throw(ArgumentError("graded channel state failed coefficient reconstruction"))
  return nothing
end

function accept_native_step!(channels, indices, step, order)
  store_corrections!(channels, indices, step.correction, order)
  store_births!(channels, step.solution.newborn, step.births, order)
  return channels
end

function native_recurrence(
  algorithm::Union{BlochFeshbach,HoriDeprit},
  representation::NativeRepresentation,
  inverse::HomologicalInverse,
  L::Harmonics,
  leading::AbstractMatrix{<:Number},
  N::Int,
  tolerance::Float64,
)
  weights = ones(eltype(leading), size(leading, 2))
  return native_recurrence(
    algorithm, representation, inverse, L, leading, weights, N, tolerance
  )
end

function native_recurrence(
  ::BlochFeshbach,
  representation::NativeRepresentation,
  inverse::HomologicalInverse,
  L::Harmonics{X},
  leading::AbstractMatrix{T},
  weights::AbstractVector{<:Number},
  N::Int,
  tolerance::Float64,
) where {X,T<:Number}
  N >= 0 || throw(ArgumentError("retained order must be nonnegative"))
  haskey(L, 0) || throw(ArgumentError("the averaged generator harmonic L_0 is required"))
  L0 = L[0]
  identity = one(L0)
  zero_element = zero(L0)
  dimension = size(leading, 1)
  channels = initial_native_channels(leading, weights)
  check_native_coefficient(representation, L0, channels, 0, dimension, tolerance)

  E = X[L0]
  H = [native_hamiltonian(representation, L0)]
  S = X[]
  slots = X[]
  known_products = Matrix{T}[]
  Y = Harmonics{X}[Harmonics{X}(0 => identity)]
  N == 0 && return NativeRecurrence(E, S, slots, H, Y, channels, known_products)
  push!(Y, harmonic_integral(L))

  for order in 1:N
    offset = intrinsic_offset(bf_truncated_kick(Y, order), S, order, identity, zero_element)
    base = -offset
    trial = harmonic_with_static(Y[order + 1], base)
    Vhat = harmonic_average(bf_static_residual(L, Y, E, order, trial), zero_element)

    known = known_gram(channels, order, dimension)
    indices, active = active_channels(channels, order, dimension)
    rates = active_weights(channels, indices)
    step = native_static_step(
      representation, inverse, L0, Vhat, known, active, rates, tolerance
    )

    push!(S, step.S)
    push!(slots, base + step.S)
    push!(E, step.E)
    push!(H, step.H)
    push!(known_products, known)

    Y[order + 1] = harmonic_with_static(Y[order + 1], slots[end])
    accept_native_step!(channels, indices, step, order)
    check_native_coefficient(representation, step.E, channels, order, dimension, tolerance)

    if order < N
      push!(Y, harmonic_integral(bf_static_residual(L, Y, E, order, Y[order + 1])))
    end
  end
  return NativeRecurrence(E, S, slots, H, Y, channels, known_products)
end

function native_recurrence(
  ::HoriDeprit,
  representation::NativeRepresentation,
  inverse::HomologicalInverse,
  L::Harmonics{X},
  leading::AbstractMatrix{T},
  weights::AbstractVector{<:Number},
  N::Int,
  tolerance::Float64,
) where {X,T<:Number}
  N >= 0 || throw(ArgumentError("retained order must be nonnegative"))
  haskey(L, 0) || throw(ArgumentError("the averaged generator harmonic L_0 is required"))
  L0 = L[0]
  identity = one(L0)
  zero_element = zero(L0)
  dimension = size(leading, 1)
  channels = initial_native_channels(leading, weights)
  check_native_coefficient(representation, L0, channels, 0, dimension, tolerance)

  E = X[L0]
  H = [native_hamiltonian(representation, L0)]
  S = X[]
  slots = X[]
  known_products = Matrix{T}[]
  G = Harmonics{X}[Harmonics{X}(), harmonic_integral(L)]
  N == 0 && return NativeRecurrence(E, S, slots, H, G, channels, known_products)

  for order in 1:N
    truncated = hd_truncated_generator(G, order)
    kick = kick_exp(truncated, identity, order)
    base = -intrinsic_offset(kick, S, order, identity, zero_element)
    G[order + 1] = harmonic_with_static(G[order + 1], base)

    Vhat = harmonic_average(hd_generator_series(L, G, order)[order + 1], zero_element)
    known = known_gram(channels, order, dimension)
    indices, active = active_channels(channels, order, dimension)
    rates = active_weights(channels, indices)
    step = native_static_step(
      representation, inverse, L0, Vhat, known, active, rates, tolerance
    )

    push!(S, step.S)
    push!(slots, base + step.S)
    push!(E, step.E)
    push!(H, step.H)
    push!(known_products, known)

    G[order + 1] = harmonic_with_static(G[order + 1], slots[end])
    accept_native_step!(channels, indices, step, order)
    check_native_coefficient(representation, step.E, channels, order, dimension, tolerance)

    W = hd_generator_series(L, G, order)[order + 1]
    native_matches(
      representation, harmonic_average(W, zero_element), step.E, 20 * tolerance
    ) || throw(
      ArgumentError("HD accepted static slot does not reproduce the retained coefficient")
    )
    order < N && push!(G, harmonic_integral(W))
  end
  return NativeRecurrence(E, S, slots, H, G, channels, known_products)
end
