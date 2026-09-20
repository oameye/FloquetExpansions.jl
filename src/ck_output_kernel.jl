struct CKOutputBlock
  vertex_start::Int
  vertex_stop::Int
  output_start::Int
  output_stop::Int
end

function Base.:(==)(left::CKOutputBlock, right::CKOutputBlock)
  return left.vertex_start == right.vertex_start &&
         left.vertex_stop == right.vertex_stop &&
         left.output_start == right.output_start &&
         left.output_stop == right.output_stop
end

Base.isequal(left::CKOutputBlock, right::CKOutputBlock) = left == right

function Base.hash(block::CKOutputBlock, seed::UInt)
  result = hash(block.vertex_start, seed)
  result = hash(block.vertex_stop, result)
  result = hash(block.output_start, result)
  return hash(block.output_stop, result)
end

struct CKOutputConstraint{H}
  block::Int
  output_stop::Int
  harmonic::H
end

function Base.:(==)(left::CKOutputConstraint, right::CKOutputConstraint)
  return left.block == right.block &&
         left.output_stop == right.output_stop &&
         isequal(left.harmonic, right.harmonic)
end

Base.isequal(left::CKOutputConstraint, right::CKOutputConstraint) = left == right

function Base.hash(constraint::CKOutputConstraint, seed::UInt)
  result = hash(constraint.block, seed)
  result = hash(constraint.output_stop, result)
  return hash(constraint.harmonic, result)
end

struct CKOutputKernelKey{H}
  output_channels::Vector{Int}
  blocks::Vector{CKOutputBlock}
  model_constraints::Vector{CKOutputConstraint{H}}
  resolvent_constraints::Vector{CKOutputConstraint{H}}
end

function Base.:(==)(left::CKOutputKernelKey, right::CKOutputKernelKey)
  return isequal(left.output_channels, right.output_channels) &&
         isequal(left.blocks, right.blocks) &&
         isequal(left.model_constraints, right.model_constraints) &&
         isequal(left.resolvent_constraints, right.resolvent_constraints)
end

Base.isequal(left::CKOutputKernelKey, right::CKOutputKernelKey) = left == right

function Base.hash(key::CKOutputKernelKey, seed::UInt)
  result = hash(length(key.output_channels), seed)
  for channel in key.output_channels
    result = hash(channel, result)
  end
  result = hash(length(key.blocks), result)
  for block in key.blocks
    result = hash(block, result)
  end
  result = hash(length(key.model_constraints), result)
  for constraint in key.model_constraints
    result = hash(constraint, result)
  end
  result = hash(length(key.resolvent_constraints), result)
  for constraint in key.resolvent_constraints
    result = hash(constraint, result)
  end
  return result
end

struct CKOutputKernel{H,T}
  terms::Dict{CKOutputKernelKey{H},T}
  zero_component::T
end

struct CKOutputPeriodPolynomial{H,T}
  coefficients::Vector{CKOutputKernel{H,T}}
end

function ck_output_prefix_data(key::CKEndpointKey{H}) where {H}
  vertex_count = length(key.vertices)
  harmonic_prefix = Vector{H}(undef, vertex_count + 1)
  harmonic_prefix[1] = zero(H)
  output_prefix = zeros(Int, vertex_count + 1)
  output_channels = Int[]

  for (vertex_index, vertex) in enumerate(key.vertices)
    harmonic_prefix[vertex_index + 1] = harmonic_prefix[vertex_index] + vertex.harmonic
    output_prefix[vertex_index + 1] = output_prefix[vertex_index]
    if ck_is_jump(vertex)
      output_prefix[vertex_index + 1] += 1
      push!(output_channels, vertex.output_channel)
    end
  end
  return (; harmonic_prefix, output_prefix, output_channels)
end

function ck_output_interval_starts(key::CKEndpointKey)
  starts = Int[]
  for interval in key.model_intervals
    interval.start in starts || push!(starts, interval.start)
  end
  for interval in key.resolvent_intervals
    interval.start in starts || push!(starts, interval.start)
  end
  sort!(starts)
  return starts
end

function ck_output_block_stop(key::CKEndpointKey, start::Int)
  stop = start
  for interval in key.model_intervals
    interval.start == start && (stop = max(stop, interval.stop))
  end
  for interval in key.resolvent_intervals
    interval.start == start && (stop = max(stop, interval.stop))
  end
  return stop
end

function ck_output_blocks(key::CKEndpointKey, output_prefix::Vector{Int})
  starts = ck_output_interval_starts(key)
  blocks = CKOutputBlock[]
  start_to_block = Dict{Int,Int}()
  previous_stop = 0

  for (block_index, start) in enumerate(starts)
    stop = ck_output_block_stop(key, start)
    block_index > 1 &&
      start < previous_stop &&
      throw(ArgumentError("endpoint constraints are not block-local prefix chains"))
    push!(
      blocks, CKOutputBlock(start, stop, output_prefix[start + 1], output_prefix[stop + 1])
    )
    start_to_block[start] = block_index
    previous_stop = stop
  end
  return blocks, start_to_block
end

function ck_output_constraint(
  interval::CKEndpointConstraint,
  block::Int,
  output_prefix::Vector{Int},
  harmonic_prefix::Vector{H},
) where {H}
  harmonic = harmonic_prefix[interval.stop + 1] - harmonic_prefix[interval.start + 1]
  return CKOutputConstraint(block, output_prefix[interval.stop + 1], harmonic)
end

function ck_output_kernel_key(key::CKEndpointKey{H}) where {H}
  prefix = ck_output_prefix_data(key)
  blocks, start_to_block = ck_output_blocks(key, prefix.output_prefix)
  model_constraints = CKOutputConstraint{H}[
    ck_output_constraint(
      interval, start_to_block[interval.start], prefix.output_prefix, prefix.harmonic_prefix
    ) for interval in key.model_intervals
  ]
  resolvent_constraints = CKOutputConstraint{H}[
    ck_output_constraint(
      interval, start_to_block[interval.start], prefix.output_prefix, prefix.harmonic_prefix
    ) for interval in key.resolvent_intervals
  ]
  return CKOutputKernelKey(
    prefix.output_channels, blocks, model_constraints, resolvent_constraints
  )
end

function ck_output_accumulate!(
  terms::Dict{CKOutputKernelKey{H},T},
  key::CKOutputKernelKey{H},
  value::T,
  zero_component::T,
) where {H,T}
  updated = get(terms, key, zero_component) + value
  if ck_kernel_iszero(updated)
    haskey(terms, key) && delete!(terms, key)
  else
    terms[key] = updated
  end
  return terms
end

function ck_output_kernel(state::CKEndpointKernel{H,T}) where {H,T}
  terms = Dict{CKOutputKernelKey{H},T}()
  for (key, value) in state.terms
    ck_output_accumulate!(terms, ck_output_kernel_key(key), value, state.zero_component)
  end
  return CKOutputKernel(terms, state.zero_component)
end

function ck_output_kernel(polynomial::CKPeriodPolynomial{H,T}) where {H,T}
  return CKOutputPeriodPolynomial(
    CKOutputKernel{H,T}[
      ck_output_kernel(coefficient) for coefficient in polynomial.coefficients
    ],
  )
end

function ck_output_constraint_mismatch(
  key::CKOutputKernelKey{H}, constraint::CKOutputConstraint{H}, sidebands::AbstractVector{H}
) where {H}
  block = key.blocks[constraint.block]
  mismatch = constraint.harmonic
  for output_index in (block.output_start + 1):constraint.output_stop
    mismatch -= sidebands[output_index]
  end
  return mismatch
end

function ck_output_ordered_sideband_coefficient(
  state::CKOutputKernel{H,T},
  output_channels::AbstractVector{Int},
  sidebands::AbstractVector{H};
  inverse_weight,
) where {H,T}
  length(output_channels) == length(sidebands) ||
    throw(ArgumentError("output channels and sidebands must have equal length"))

  result = state.zero_component
  for (key, value) in state.terms
    key.output_channels == output_channels || continue

    valid = true
    for constraint in key.model_constraints
      if !iszero(ck_output_constraint_mismatch(key, constraint, sidebands))
        valid = false
        break
      end
    end
    valid || continue

    weighted_value = value
    for constraint in key.resolvent_constraints
      mismatch = ck_output_constraint_mismatch(key, constraint, sidebands)
      if iszero(mismatch)
        valid = false
        break
      end
      weighted_value = inverse_weight(mismatch) * weighted_value
    end
    valid || continue
    result += weighted_value
  end
  return result
end

function ck_output_ordered_sideband_coefficients(
  polynomial::CKOutputPeriodPolynomial{H,T},
  output_channels::AbstractVector{Int},
  sidebands::AbstractVector{H};
  inverse_weight,
) where {H,T}
  return T[
    ck_output_ordered_sideband_coefficient(
      coefficient, output_channels, sidebands; inverse_weight=inverse_weight
    ) for coefficient in polynomial.coefficients
  ]
end

function ck_output_kernel_complete(key::CKOutputKernelKey)
  isempty(key.output_channels) && return true
  constrained = falses(length(key.output_channels))
  for constraint in key.model_constraints
    iszero(constraint.output_stop) || (constrained[constraint.output_stop] = true)
  end
  for constraint in key.resolvent_constraints
    iszero(constraint.output_stop) || (constrained[constraint.output_stop] = true)
  end
  return all(constrained)
end
