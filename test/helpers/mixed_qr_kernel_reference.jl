const MixedKernelExact = Complex{Rational{Int}}
const mixed_kernel_im = MixedKernelExact(0 // 1, 1 // 1)

struct MixedQRVertex
  kind::Symbol
  channel::Symbol
  harmonic::Int
end

mixed_qr_jump(channel::Symbol, harmonic::Int) = MixedQRVertex(:jump, channel, harmonic)
mixed_qr_drift(harmonic::Int) = MixedQRVertex(:drift, :drift, harmonic)
mixed_qr_is_jump(vertex::MixedQRVertex) = vertex.kind === :jump

function mixed_add_coefficient!(
  state::Dict{Int,Vector{MixedKernelExact}},
  frequency::Int,
  degree::Int,
  value::MixedKernelExact,
)
  coefficients = get!(state, frequency) do
    return MixedKernelExact[]
  end
  while length(coefficients) <= degree
    push!(coefficients, zero(MixedKernelExact))
  end
  coefficients[degree + 1] += value
  return state
end

function mixed_scale_add!(
  target::Dict{Int,Vector{MixedKernelExact}},
  source::Dict{Int,Vector{MixedKernelExact}},
  scale::MixedKernelExact,
)
  for (frequency, coefficients) in source
    for (degree_index, coefficient) in enumerate(coefficients)
      mixed_add_coefficient!(target, frequency, degree_index - 1, scale * coefficient)
    end
  end
  return target
end

function mixed_monomial_integral(frequency::Int, degree::Int)
  degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
  result = Dict{Int,Vector{MixedKernelExact}}()

  if iszero(frequency)
    coefficient = MixedKernelExact(1 // (degree + 1), 0 // 1)
    mixed_add_coefficient!(result, 0, degree + 1, coefficient)
    return result
  end

  inverse_frequency = one(MixedKernelExact) / (mixed_kernel_im * frequency)
  if iszero(degree)
    mixed_add_coefficient!(result, 0, 0, inverse_frequency)
    mixed_add_coefficient!(result, frequency, 0, -inverse_frequency)
    return result
  end

  lower = mixed_monomial_integral(frequency, degree - 1)
  scale = MixedKernelExact(degree // 1, 0 // 1) * inverse_frequency
  mixed_scale_add!(result, lower, scale)
  mixed_add_coefficient!(result, frequency, degree, -inverse_frequency)
  return result
end

function mixed_shift_frequency(state::Dict{Int,Vector{MixedKernelExact}}, harmonic::Int)
  shifted = Dict{Int,Vector{MixedKernelExact}}()
  for (frequency, coefficients) in state
    shifted[frequency + harmonic] = copy(coefficients)
  end
  return shifted
end

function mixed_integrate_fourier_state(state::Dict{Int,Vector{MixedKernelExact}})
  integrated = Dict{Int,Vector{MixedKernelExact}}()
  for (frequency, coefficients) in state
    for (degree_index, coefficient) in enumerate(coefficients)
      primitive = mixed_monomial_integral(frequency, degree_index - 1)
      mixed_scale_add!(integrated, primitive, coefficient)
    end
  end
  return integrated
end

function mixed_trim_polynomial!(polynomial::Vector{MixedKernelExact})
  while length(polynomial) > 1 && iszero(last(polynomial))
    pop!(polynomial)
  end
  return polynomial
end

function mixed_add_polynomials(
  left::Vector{MixedKernelExact}, right::Vector{MixedKernelExact}
)
  degree = max(length(left), length(right))
  result = fill(zero(MixedKernelExact), degree)
  for index in eachindex(left)
    result[index] += left[index]
  end
  for index in eachindex(right)
    result[index] += right[index]
  end
  return mixed_trim_polynomial!(result)
end

function mixed_formal_period_simplex(frequencies)
  state = Dict{Int,Vector{MixedKernelExact}}(0 => MixedKernelExact[one(MixedKernelExact)])
  for frequency in frequencies
    state = mixed_integrate_fourier_state(mixed_shift_frequency(state, frequency))
  end

  polynomial = MixedKernelExact[zero(MixedKernelExact)]
  for coefficients in values(state)
    for (degree_index, coefficient) in enumerate(coefficients)
      while length(polynomial) < degree_index
        push!(polynomial, zero(MixedKernelExact))
      end
      polynomial[degree_index] += coefficient
    end
  end
  return mixed_trim_polynomial!(polynomial)
end

function mixed_period_coefficient(polynomial::Vector{MixedKernelExact}, degree::Int)
  degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
  degree + 1 <= length(polynomial) || return zero(MixedKernelExact)
  return polynomial[degree + 1]
end

function mixed_qr_jump_vertices(word::Vector{MixedQRVertex})
  return MixedQRVertex[vertex for vertex in word if mixed_qr_is_jump(vertex)]
end

function mixed_qr_add_precedence!(
  predecessors::Vector{Vector{Int}}, before::Int, after::Int
)
  before == after && return predecessors
  before in predecessors[after] || push!(predecessors[after], before)
  return predecessors
end

function mixed_qr_poset(left::Vector{MixedQRVertex}, right::Vector{MixedQRVertex})
  left_jumps = mixed_qr_jump_vertices(left)
  right_jumps = mixed_qr_jump_vertices(right)
  length(left_jumps) == length(right_jumps) || return nothing
  all(
    left_jump.channel == right_jump.channel for
    (left_jump, right_jump) in zip(left_jumps, right_jumps)
  ) || return nothing

  output_number = length(left_jumps)
  frequencies = Int[
    left_jumps[index].harmonic - right_jumps[index].harmonic for index in 1:output_number
  ]
  predecessors = [Int[] for _ in 1:output_number]

  left_nodes = Int[]
  jump_index = 0
  for vertex in left
    node = if mixed_qr_is_jump(vertex)
      jump_index += 1
      jump_index
    else
      push!(frequencies, vertex.harmonic)
      push!(predecessors, Int[])
      length(frequencies)
    end
    push!(left_nodes, node)
  end

  right_nodes = Int[]
  jump_index = 0
  for vertex in right
    node = if mixed_qr_is_jump(vertex)
      jump_index += 1
      jump_index
    else
      push!(frequencies, -vertex.harmonic)
      push!(predecessors, Int[])
      length(frequencies)
    end
    push!(right_nodes, node)
  end

  for nodes in (left_nodes, right_nodes)
    for index in 1:(length(nodes) - 1)
      mixed_qr_add_precedence!(predecessors, nodes[index], nodes[index + 1])
    end
  end

  return (; frequencies, predecessors, output_number)
end

function mixed_qr_linear_extensions(predecessors::Vector{Vector{Int}})
  node_count = length(predecessors)
  used = falses(node_count)
  current = Int[]
  extensions = Vector{Vector{Int}}()

  function visit!()
    if length(current) == node_count
      push!(extensions, copy(current))
      return nothing
    end

    for node in 1:node_count
      used[node] && continue
      all(predecessor -> used[predecessor], predecessors[node]) || continue
      used[node] = true
      push!(current, node)
      visit!()
      pop!(current)
      used[node] = false
    end
  end

  visit!()
  return extensions
end

function mixed_qr_gram_polynomial(left::Vector{MixedQRVertex}, right::Vector{MixedQRVertex})
  poset = mixed_qr_poset(left, right)
  isnothing(poset) && return MixedKernelExact[zero(MixedKernelExact)]

  extensions = mixed_qr_linear_extensions(poset.predecessors)
  isempty(extensions) && error("mixed Q/R branch constraints contain a cycle")

  result = MixedKernelExact[zero(MixedKernelExact)]
  for extension in extensions
    frequencies = Int[poset.frequencies[node] for node in extension]
    result = mixed_add_polynomials(result, mixed_formal_period_simplex(frequencies))
  end
  return result
end

function mixed_one_drift_jump_kernel(
  position::Symbol, drift_harmonic::Int, jump_harmonic::Int
)
  iszero(drift_harmonic) && throw(ArgumentError("drift harmonic must be nonzero"))
  position in (:before, :after) ||
    throw(ArgumentError("position must be :before or :after"))

  inverse_frequency = one(MixedKernelExact) / (mixed_kernel_im * drift_harmonic)
  leading = position === :before ? inverse_frequency : -inverse_frequency
  return Dict(jump_harmonic => leading, jump_harmonic + drift_harmonic => -leading)
end

function mixed_qr_word_count(jump_labels::Int, drift_labels::Int, delta_cutoff::Int)
  jump_labels >= 1 || throw(ArgumentError("jump label count must be positive"))
  drift_labels >= 1 || throw(ArgumentError("drift label count must be positive"))
  delta_cutoff >= 0 || throw(ArgumentError("delta cutoff must be nonnegative"))

  total = 0
  for jump_count in 0:delta_cutoff
    max_drifts = div(delta_cutoff - jump_count, 2)
    for drift_count in 0:max_drifts
      vertex_count = jump_count + drift_count
      total +=
        binomial(vertex_count, jump_count) *
        jump_labels^jump_count *
        drift_labels^drift_count
    end
  end
  return total
end
