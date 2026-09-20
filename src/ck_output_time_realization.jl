struct CKOutputTimeCoordinate{H,T}
  block::Int
  output_stop::Int
  kernel::CKResolventTimeKernel{H,T}
end

struct CKOutputTimeRealization{H,T}
  blocks::Vector{CKOutputBlock}
  coordinates::Vector{CKOutputTimeCoordinate{H,T}}
  scalar::T
  valid::Bool
end

function ck_phase_time_kernel(harmonic::H, coefficient::T) where {H<:Integer,T}
  return CKResolventTimeKernel(
    CKHomologicalComponent{H,T}[],
    CKPhaseComponent{H,T}[CKPhaseComponent(harmonic, coefficient)],
    zero(coefficient),
  )
end

function ck_output_constraint_harmonics(
  constraints::Vector{CKOutputConstraint{H}}, block::Int, output_stop::Int
) where {H}
  result = H[]
  for constraint in constraints
    constraint.block == block || continue
    constraint.output_stop == output_stop || continue
    push!(result, constraint.harmonic)
  end
  return result
end

function ck_output_unique_model_harmonic(model_harmonics::Vector{H}) where {H<:Integer}
  isempty(model_harmonics) && return zero(H), false, true
  harmonic = first(model_harmonics)
  consistent = all(candidate -> isequal(candidate, harmonic), model_harmonics)
  return harmonic, true, consistent
end

function ck_output_scalar_constraint_weight(
  model_harmonics::Vector{H}, resolvent_harmonics::Vector{H}, imaginary::T
) where {H<:Integer,T}
  zero_coefficient = zero(imaginary)
  model_harmonic, model_present, model_consistent =
    ck_output_unique_model_harmonic(model_harmonics)
  model_consistent || return zero_coefficient, false
  model_present && !iszero(model_harmonic) && return zero_coefficient, false

  weight = one(imaginary)
  for harmonic in resolvent_harmonics
    iszero(harmonic) && return zero_coefficient, false
    weight *= imaginary / harmonic
  end
  return weight, true
end

function ck_output_variable_constraint_kernel(
  model_harmonics::Vector{H}, resolvent_harmonics::Vector{H}, imaginary::T
) where {H<:Integer,T}
  zero_coefficient = zero(imaginary)
  one_coefficient = one(imaginary)
  model_harmonic, model_present, model_consistent =
    ck_output_unique_model_harmonic(model_harmonics)
  model_consistent ||
    return ck_phase_time_kernel(zero(H), zero_coefficient), zero_coefficient, false

  if model_present
    scalar = one_coefficient
    for harmonic in resolvent_harmonics
      mismatch = harmonic - model_harmonic
      iszero(mismatch) &&
        return ck_phase_time_kernel(model_harmonic, zero_coefficient), zero_coefficient, false
      scalar *= imaginary / mismatch
    end
    return ck_phase_time_kernel(model_harmonic, one_coefficient), scalar, true
  end

  isempty(resolvent_harmonics) &&
    return ck_phase_time_kernel(zero(H), zero_coefficient), zero_coefficient, false
  return ck_resolvent_time_kernel(resolvent_harmonics, imaginary), one_coefficient, true
end

function ck_output_time_realization(
  key::CKOutputKernelKey{H}, imaginary::T
) where {H<:Integer,T}
  scalar = one(imaginary)
  coordinates = CKOutputTimeCoordinate{H,T}[]

  for (block_index, block) in enumerate(key.blocks)
    scalar_models = ck_output_constraint_harmonics(
      key.model_constraints, block_index, block.output_start
    )
    scalar_resolvents = ck_output_constraint_harmonics(
      key.resolvent_constraints, block_index, block.output_start
    )
    scalar_weight, scalar_valid = ck_output_scalar_constraint_weight(
      scalar_models, scalar_resolvents, imaginary
    )
    scalar_valid ||
      return CKOutputTimeRealization(key.blocks, coordinates, zero(imaginary), false)
    scalar *= scalar_weight

    for output_stop in (block.output_start + 1):block.output_stop
      model_harmonics = ck_output_constraint_harmonics(
        key.model_constraints, block_index, output_stop
      )
      resolvent_harmonics = ck_output_constraint_harmonics(
        key.resolvent_constraints, block_index, output_stop
      )
      coordinate_kernel, coordinate_scalar, coordinate_valid =
        ck_output_variable_constraint_kernel(model_harmonics, resolvent_harmonics, imaginary)
      coordinate_valid ||
        return CKOutputTimeRealization(key.blocks, coordinates, zero(imaginary), false)
      scalar *= coordinate_scalar
      push!(coordinates, CKOutputTimeCoordinate(block_index, output_stop, coordinate_kernel))
    end
  end

  length(coordinates) == length(key.output_channels) ||
    return CKOutputTimeRealization(key.blocks, coordinates, zero(imaginary), false)
  return CKOutputTimeRealization(key.blocks, coordinates, scalar, true)
end

function ck_output_time_cumulative_sideband(
  realization::CKOutputTimeRealization{H},
  coordinate::CKOutputTimeCoordinate{H},
  sidebands::AbstractVector{H},
) where {H<:Integer}
  block = realization.blocks[coordinate.block]
  result = zero(H)
  for output_index in (block.output_start + 1):coordinate.output_stop
    result += sidebands[output_index]
  end
  return result
end

function ck_output_time_fourier_weight(
  realization::CKOutputTimeRealization{H,T}, sidebands::AbstractVector{H}, imaginary::T
) where {H<:Integer,T}
  realization.valid || return zero(imaginary)
  result = realization.scalar
  for coordinate in realization.coordinates
    cumulative_sideband = ck_output_time_cumulative_sideband(
      realization, coordinate, sidebands
    )
    result *= ck_resolvent_time_fourier_coefficient(
      coordinate.kernel, cumulative_sideband, imaginary
    )
  end
  return result
end
