ck_channel_gauge_commutator(left, right) = left * right - right * left

function ck_channel_gauge_hamiltonian_action(H)
  identity_matrix = Matrix{eltype(H)}(I, size(H, 1), size(H, 2))
  return -ck_channel_gauge_im * (
    kron(identity_matrix, H) - kron(transpose(H), identity_matrix)
  )
end

function ck_channel_gauge_cross_dissipator(left, right)
  identity_matrix = Matrix{eltype(left)}(I, size(left, 1), size(left, 2))
  norm = adjoint(right) * left
  return kron(conj(right), left) -
         (1 // 2) * kron(identity_matrix, norm) -
         (1 // 2) * kron(transpose(norm), identity_matrix)
end

function ck_channel_gauge_accumulate!(terms, harmonic, value, zero_component)
  updated = get(terms, harmonic, zero_component) + value
  if all(iszero, updated)
    haskey(terms, harmonic) && delete!(terms, harmonic)
  else
    terms[harmonic] = updated
  end
  return terms
end

function ck_channel_gauge_loss_harmonics(jumps, zero_component)
  result = Dict{Int,typeof(zero_component)}()
  for (left_harmonic, left) in jumps, (right_harmonic, right) in jumps
    harmonic = left_harmonic - right_harmonic
    contribution = -(1 // 2) * adjoint(right) * left
    ck_channel_gauge_accumulate!(result, harmonic, contribution, zero_component)
  end
  return result
end

function ck_channel_gauge_dissipator_harmonics(jumps, zero_component)
  zero_superoperator = zeros(eltype(zero_component), size(zero_component, 1)^2, size(zero_component, 2)^2)
  result = Dict{Int,typeof(zero_superoperator)}()
  for (left_harmonic, left) in jumps, (right_harmonic, right) in jumps
    harmonic = left_harmonic - right_harmonic
    contribution = ck_channel_gauge_cross_dissipator(left, right)
    ck_channel_gauge_accumulate!(result, harmonic, contribution, zero_superoperator)
  end
  return result
end

function ck_channel_gauge_first_kicks(hamiltonian)
  result = Dict{Int,Matrix{CKChannelGaugeExact}}()
  for (harmonic, value) in hamiltonian
    iszero(harmonic) && continue
    result[harmonic] = ck_channel_gauge_im * (1 // harmonic) * value
  end
  return result
end

function ck_channel_gauge_second_kicks(hamiltonian, zero_component)
  result = Dict{Int,Matrix{CKChannelGaugeExact}}()
  H0 = hamiltonian[0]
  nonzero_harmonics = filter(!iszero, collect(keys(hamiltonian)))
  for harmonic in nonzero_harmonics
    value =
      -ck_channel_gauge_im * (1 // harmonic^2) *
      ck_channel_gauge_commutator(hamiltonian[harmonic], H0)
    for inner_harmonic in nonzero_harmonics
      inner_harmonic == harmonic && continue
      nested = get(hamiltonian, harmonic - inner_harmonic, zero_component)
      all(iszero, nested) && continue
      value +=
        -ck_channel_gauge_im * (1 // (2 * harmonic * inner_harmonic)) *
        ck_channel_gauge_commutator(hamiltonian[inner_harmonic], nested)
    end
    all(iszero, value) || (result[harmonic] = value)
  end
  return result
end

function ck_channel_gauge_transported_jump_orders(hamiltonian, jumps, zero_component)
  first_kicks = ck_channel_gauge_first_kicks(hamiltonian)
  second_kicks = ck_channel_gauge_second_kicks(hamiltonian, zero_component)
  first = Dict{Int,Matrix{CKChannelGaugeExact}}()
  second = Dict{Int,Matrix{CKChannelGaugeExact}}()

  for (kick_harmonic, kick) in first_kicks, (jump_harmonic, jump) in jumps
    harmonic = kick_harmonic + jump_harmonic
    contribution = ck_channel_gauge_im * ck_channel_gauge_commutator(kick, jump)
    ck_channel_gauge_accumulate!(first, harmonic, contribution, zero_component)
  end

  for (kick_harmonic, kick) in second_kicks, (jump_harmonic, jump) in jumps
    harmonic = kick_harmonic + jump_harmonic
    contribution = ck_channel_gauge_im * ck_channel_gauge_commutator(kick, jump)
    ck_channel_gauge_accumulate!(second, harmonic, contribution, zero_component)
  end
  for (left_harmonic, left_kick) in first_kicks,
    (right_harmonic, right_kick) in first_kicks,
    (jump_harmonic, jump) in jumps
    harmonic = left_harmonic + right_harmonic + jump_harmonic
    contribution =
      -(1 // 2) *
      ck_channel_gauge_commutator(
        left_kick, ck_channel_gauge_commutator(right_kick, jump)
      )
    ck_channel_gauge_accumulate!(second, harmonic, contribution, zero_component)
  end
  return first, second
end

function ck_channel_gauge_cp_second_order(hamiltonian, jumps, zero_component)
  first, second = ck_channel_gauge_transported_jump_orders(
    hamiltonian, jumps, zero_component
  )
  zero_superoperator = zeros(CKChannelGaugeExact, 4, 4)
  result = copy(zero_superoperator)
  sidebands = union(keys(jumps), keys(first), keys(second))
  for harmonic in sidebands
    bare = get(jumps, harmonic, zero_component)
    first_order = get(first, harmonic, zero_component)
    second_order = get(second, harmonic, zero_component)
    result += ck_channel_gauge_cross_dissipator(bare, second_order)
    result += ck_channel_gauge_cross_dissipator(second_order, bare)
    result += ck_channel_gauge_cross_dissipator(first_order, first_order)
  end
  return result
end

function ck_channel_gauge_static_correction(hamiltonian, jumps, zero_component)
  H_map = Dict(
    harmonic => ck_channel_gauge_hamiltonian_action(value) for
    (harmonic, value) in hamiltonian
  )
  R = ck_channel_gauge_dissipator_harmonics(jumps, zero_component)
  zero_superoperator = zeros(CKChannelGaugeExact, 4, 4)
  correction = copy(zero_superoperator)
  for harmonic in keys(H_map)
    iszero(harmonic) && continue
    correction +=
      (1 // (2 * harmonic^2)) *
      ck_channel_gauge_commutator(
        H_map[harmonic], get(R, -harmonic, zero_superoperator)
      )
  end
  return correction, H_map[0]
end
