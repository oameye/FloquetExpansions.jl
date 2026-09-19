@enum CKKernelSector::UInt8 begin
  CKUnresolvedSector = 0x01
  CKModelSector = 0x02
  CKComplementSector = 0x03
end

struct CKBranchVertex{H}
  output_channel::Int
  harmonic::H
end

function ck_jump_vertex(output_channel::Int, harmonic::H) where {H}
  output_channel >= 1 || throw(ArgumentError("output channel must be positive"))
  return CKBranchVertex(output_channel, harmonic)
end

ck_drift_vertex(harmonic::H) where {H} = CKBranchVertex(0, harmonic)
ck_is_jump(vertex::CKBranchVertex) = vertex.output_channel > 0

struct CKKernelKey{H}
  sector::CKKernelSector
  vertices::Vector{CKBranchVertex{H}}
  model_cuts::Vector{Int}
  resolvent_cuts::Vector{Int}
end

function CKKernelKey(
  sector::CKKernelSector,
  vertices::Vector{CKBranchVertex{H}},
  model_cuts::Vector{Int},
  resolvent_cuts::Vector{Int},
) where {H}
  vertex_count = length(vertices)
  all(cut -> 1 <= cut <= vertex_count, model_cuts) ||
    throw(ArgumentError("model cuts must lie inside the physical branch"))
  all(cut -> 1 <= cut <= vertex_count, resolvent_cuts) ||
    throw(ArgumentError("resolvent cuts must lie inside the physical branch"))

  canonical_model_cuts = unique(sort(model_cuts))
  canonical_resolvent_cuts = sort(resolvent_cuts)
  return CKKernelKey{H}(
    sector, copy(vertices), canonical_model_cuts, canonical_resolvent_cuts
  )
end

function Base.:(==)(left::CKKernelKey, right::CKKernelKey)
  return left.sector == right.sector &&
         left.vertices == right.vertices &&
         left.model_cuts == right.model_cuts &&
         left.resolvent_cuts == right.resolvent_cuts
end

Base.isequal(left::CKKernelKey, right::CKKernelKey) = left == right

function Base.hash(key::CKKernelKey, seed::UInt)
  result = hash(key.sector, seed)
  result = hash(length(key.vertices), result)
  for vertex in key.vertices
    result = hash(vertex, result)
  end
  result = hash(length(key.model_cuts), result)
  for cut in key.model_cuts
    result = hash(cut, result)
  end
  result = hash(length(key.resolvent_cuts), result)
  for cut in key.resolvent_cuts
    result = hash(cut, result)
  end
  return result
end

struct CKPhysicalKernel{H,T}
  terms::Dict{CKKernelKey{H},T}
  zero_component::T
end

ck_kernel_iszero(value) = iszero(value)
ck_kernel_iszero(value::AbstractArray) = all(iszero, value)

function ck_kernel_output_number(key::CKKernelKey)
  return count(ck_is_jump, key.vertices)
end

function ck_kernel_drift_number(key::CKKernelKey)
  return length(key.vertices) - ck_kernel_output_number(key)
end

function ck_kernel_total_harmonic(key::CKKernelKey{H}) where {H}
  total = zero(H)
  for vertex in key.vertices
    total += vertex.harmonic
  end
  return total
end

function ck_kernel_has_forced_zero(key::CKKernelKey{H}) where {H}
  any(cut -> cut in key.resolvent_cuts, key.model_cuts) && return true

  mismatch = zero(H)
  outputs_seen = 0
  model_index = 1
  resolvent_index = 1
  for (vertex_index, vertex) in enumerate(key.vertices)
    mismatch += vertex.harmonic
    outputs_seen += Int(ck_is_jump(vertex))

    while model_index <= length(key.model_cuts) &&
          key.model_cuts[model_index] == vertex_index
      outputs_seen == 0 && !iszero(mismatch) && return true
      model_index += 1
    end
    while resolvent_index <= length(key.resolvent_cuts) &&
          key.resolvent_cuts[resolvent_index] == vertex_index
      outputs_seen == 0 && iszero(mismatch) && return true
      resolvent_index += 1
    end
  end

  if outputs_seen == 0
    key.sector == CKModelSector && !iszero(mismatch) && return true
    key.sector == CKComplementSector && iszero(mismatch) && return true
  end
  return false
end

function ck_kernel_accumulate!(
  terms::Dict{CKKernelKey{H},T}, key::CKKernelKey{H}, value::T, zero_component::T
) where {H,T}
  ck_kernel_has_forced_zero(key) && return terms
  updated = get(terms, key, zero_component) + value
  if ck_kernel_iszero(updated)
    haskey(terms, key) && delete!(terms, key)
  else
    terms[key] = updated
  end
  return terms
end

function ck_physical_kernel(terms::Dict{CKKernelKey{H},T}, zero_component::T) where {H,T}
  result = Dict{CKKernelKey{H},T}()
  for (key, value) in terms
    ck_kernel_iszero(value) || ck_kernel_accumulate!(result, key, value, zero_component)
  end
  return CKPhysicalKernel(result, zero_component)
end

function ck_kernel_zero(zero_harmonic::H, zero_component::T) where {H,T}
  return CKPhysicalKernel(Dict{CKKernelKey{H},T}(), zero_component)
end

function ck_kernel_identity(
  zero_harmonic::H, identity_component::T, zero_component::T
) where {H,T}
  terms = Dict{CKKernelKey{H},T}()
  key = CKKernelKey(CKModelSector, CKBranchVertex{H}[], Int[], Int[])
  terms[key] = identity_component
  return CKPhysicalKernel(terms, zero_component)
end

function ck_kernel_generator(
  components::AbstractDict{CKBranchVertex{H},T}, zero_component::T
) where {H,T}
  terms = Dict{CKKernelKey{H},T}()
  for (vertex, value) in components
    key = CKKernelKey(CKUnresolvedSector, CKBranchVertex{H}[vertex], Int[], Int[])
    ck_kernel_accumulate!(terms, key, value, zero_component)
  end
  return CKPhysicalKernel(terms, zero_component)
end

function Base.:(==)(left::CKPhysicalKernel, right::CKPhysicalKernel)
  return left.zero_component == right.zero_component && left.terms == right.terms
end

function Base.:+(left::CKPhysicalKernel{H,T}, right::CKPhysicalKernel{H,T}) where {H,T}
  result = copy(left.terms)
  for (key, value) in right.terms
    ck_kernel_accumulate!(result, key, value, left.zero_component)
  end
  return CKPhysicalKernel(result, left.zero_component)
end

function Base.:-(left::CKPhysicalKernel{H,T}, right::CKPhysicalKernel{H,T}) where {H,T}
  result = copy(left.terms)
  for (key, value) in right.terms
    ck_kernel_accumulate!(result, key, -value, left.zero_component)
  end
  return CKPhysicalKernel(result, left.zero_component)
end

function ck_kernel_product_sector(left::CKKernelSector, right::CKKernelSector)
  left == CKUnresolvedSector && return CKUnresolvedSector
  left == CKModelSector && right == CKModelSector && return CKModelSector
  left == CKComplementSector && right == CKModelSector && return CKComplementSector
  return nothing
end

function ck_kernel_shifted_cuts(right_cuts::Vector{Int}, left_cuts::Vector{Int}, shift::Int)
  result = copy(right_cuts)
  append!(result, (cut + shift for cut in left_cuts))
  return result
end

function ck_kernel_product(
  left::CKPhysicalKernel{H,T}, right::CKPhysicalKernel{H,T}
) where {H,T}
  result = Dict{CKKernelKey{H},T}()
  for (left_key, left_value) in left.terms, (right_key, right_value) in right.terms
    sector = ck_kernel_product_sector(left_key.sector, right_key.sector)
    isnothing(sector) && continue

    right_length = length(right_key.vertices)
    vertices = CKBranchVertex{H}[right_key.vertices; left_key.vertices]
    model_cuts = ck_kernel_shifted_cuts(
      right_key.model_cuts, left_key.model_cuts, right_length
    )
    resolvent_cuts = ck_kernel_shifted_cuts(
      right_key.resolvent_cuts, left_key.resolvent_cuts, right_length
    )
    key = CKKernelKey(sector, vertices, model_cuts, resolvent_cuts)
    ck_kernel_accumulate!(result, key, left_value * right_value, left.zero_component)
  end
  return CKPhysicalKernel(result, left.zero_component)
end

function ck_kernel_project_model(state::CKPhysicalKernel{H,T}) where {H,T}
  result = Dict{CKKernelKey{H},T}()
  for (key, value) in state.terms
    key.sector == CKComplementSector && continue

    model_key = if key.sector == CKModelSector
      key
    else
      model_cuts = copy(key.model_cuts)
      isempty(key.vertices) || push!(model_cuts, length(key.vertices))
      CKKernelKey(CKModelSector, key.vertices, model_cuts, key.resolvent_cuts)
    end
    ck_kernel_accumulate!(result, model_key, value, state.zero_component)
  end
  return CKPhysicalKernel(result, state.zero_component)
end

function ck_kernel_solve_complement(state::CKPhysicalKernel{H,T}) where {H,T}
  result = Dict{CKKernelKey{H},T}()
  for (key, value) in state.terms
    key.sector == CKModelSector && continue

    resolvent_cuts = copy(key.resolvent_cuts)
    push!(resolvent_cuts, length(key.vertices))
    solved_key = CKKernelKey(
      CKComplementSector, key.vertices, key.model_cuts, resolvent_cuts
    )
    ck_kernel_accumulate!(result, solved_key, value, state.zero_component)
  end
  return CKPhysicalKernel(result, state.zero_component)
end

function ck_kernel_ordered_sideband_coefficient(
  state::CKPhysicalKernel{H,T}, sidebands::AbstractVector{H}; inverse_weight
) where {H,T}
  result = state.zero_component
  for (key, value) in state.terms
    ck_kernel_output_number(key) == length(sidebands) || continue

    mismatch = zero(H)
    sideband_index = 0
    model_index = 1
    resolvent_index = 1
    weighted_value = value
    valid = true

    for (vertex_index, vertex) in enumerate(key.vertices)
      mismatch += vertex.harmonic
      if ck_is_jump(vertex)
        sideband_index += 1
        mismatch -= sidebands[sideband_index]
      end

      while model_index <= length(key.model_cuts) &&
            key.model_cuts[model_index] == vertex_index
        if !iszero(mismatch)
          valid = false
          break
        end
        model_index += 1
      end
      valid || break

      while resolvent_index <= length(key.resolvent_cuts) &&
            key.resolvent_cuts[resolvent_index] == vertex_index
        if iszero(mismatch)
          valid = false
          break
        end
        weighted_value = inverse_weight(mismatch) * weighted_value
        resolvent_index += 1
      end
      valid || break
    end
    valid || continue

    key.sector == CKModelSector && !iszero(mismatch) && continue
    key.sector == CKComplementSector && iszero(mismatch) && continue
    result += weighted_value
  end
  return result
end
