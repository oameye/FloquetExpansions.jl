using Test
using LinearAlgebra: I

include(joinpath(@__DIR__, "ck_periodic_bernoulli_kernel_reference.jl"))
include(joinpath(@__DIR__, "helpers", "bloch_order_kernel_reference.jl"))

const CKPhysicalB2Exact = CKBernoulliExact
const ck_physical_b2_im = CKPhysicalB2Exact(0 // 1, 1 // 1)

@enum CKPhysicalB2Kind begin
  CKPhysicalB2Identity
  CKPhysicalB2Bare
  CKPhysicalB2Model
  CKPhysicalB2Solved
  CKPhysicalB2ProjectableTwo
  CKPhysicalB2ComplementTwo
  CKPhysicalB2PhysicalTwo
end

struct CKPhysicalB2State{T}
  terms::Dict{Tuple{CKPhysicalB2Kind,Int,Int},T}
  zero_component::T
end

ck_physical_b2_iszero(value) = iszero(value)
ck_physical_b2_iszero(value::AbstractMatrix) = all(iszero, value)

function ck_physical_b2_state(terms, zero_component::T) where {T}
  cleaned = Dict{Tuple{CKPhysicalB2Kind,Int,Int},T}()
  for (key, value) in terms
    ck_physical_b2_iszero(value) || (cleaned[key] = value)
  end
  return CKPhysicalB2State(cleaned, zero_component)
end

function Base.:+(left::CKPhysicalB2State{T}, right::CKPhysicalB2State{T}) where {T}
  result = copy(left.terms)
  for (key, value) in right.terms
    updated = get(result, key, left.zero_component) + value
    if ck_physical_b2_iszero(updated)
      haskey(result, key) && delete!(result, key)
    else
      result[key] = updated
    end
  end
  return CKPhysicalB2State(result, left.zero_component)
end

function Base.:-(left::CKPhysicalB2State{T}, right::CKPhysicalB2State{T}) where {T}
  result = copy(left.terms)
  for (key, value) in right.terms
    updated = get(result, key, left.zero_component) - value
    if ck_physical_b2_iszero(updated)
      haskey(result, key) && delete!(result, key)
    else
      result[key] = updated
    end
  end
  return CKPhysicalB2State(result, left.zero_component)
end

function ck_physical_b2_bare(amplitudes::AbstractDict{Int,T}, zero_component::T) where {T}
  terms = Dict{Tuple{CKPhysicalB2Kind,Int,Int},T}()
  for (harmonic, value) in amplitudes
    terms[(CKPhysicalB2Bare, harmonic, 0)] = value
  end
  return ck_physical_b2_state(terms, zero_component)
end

function ck_physical_b2_identity(zero_component::T, identity_component::T) where {T}
  return CKPhysicalB2State(
    Dict((CKPhysicalB2Identity, 0, 0) => identity_component), zero_component
  )
end

function ck_physical_b2_product(
  left::CKPhysicalB2State{T}, right::CKPhysicalB2State{T}
) where {T}
  result = Dict{Tuple{CKPhysicalB2Kind,Int,Int},T}()

  for ((left_kind, left_first, left_second), left_value) in left.terms,
    ((right_kind, right_first, right_second), right_value) in right.terms

    key = if right_kind == CKPhysicalB2Identity
      (left_kind, left_first, left_second)
    elseif left_kind == CKPhysicalB2Identity
      (right_kind, right_first, right_second)
    elseif left_kind == CKPhysicalB2Bare && right_kind == CKPhysicalB2Solved
      # Matrix multiplication is time ordered: the later bare vertex acts on the left.
      # Store the physical output harmonics in chronological order (earlier, later).
      (CKPhysicalB2ProjectableTwo, right_first, left_first)
    elseif left_kind == CKPhysicalB2Solved && right_kind == CKPhysicalB2Model
      # This fold retains the nonzero internal mismatch of the solved leg and is
      # therefore removed by the model-space projector at this order.
      (CKPhysicalB2ComplementTwo, left_first, right_first)
    else
      throw(
        ArgumentError(
          "unsupported physical B2 product $left_kind * $right_kind in order-two reference",
        ),
      )
    end

    updated = get(result, key, left.zero_component) + left_value * right_value
    if ck_physical_b2_iszero(updated)
      haskey(result, key) && delete!(result, key)
    else
      result[key] = updated
    end
  end

  return CKPhysicalB2State(result, left.zero_component)
end

function ck_physical_b2_project_model(state::CKPhysicalB2State{T}) where {T}
  result = Dict{Tuple{CKPhysicalB2Kind,Int,Int},T}()
  for ((kind, first_harmonic, second_harmonic), value) in state.terms
    key = if kind == CKPhysicalB2Bare
      (CKPhysicalB2Model, first_harmonic, 0)
    elseif kind == CKPhysicalB2ProjectableTwo
      (CKPhysicalB2PhysicalTwo, first_harmonic, second_harmonic)
    elseif kind == CKPhysicalB2Model || kind == CKPhysicalB2PhysicalTwo
      (kind, first_harmonic, second_harmonic)
    else
      continue
    end
    result[key] = get(result, key, state.zero_component) + value
  end
  return ck_physical_b2_state(result, state.zero_component)
end

function ck_physical_b2_solve_complement(state::CKPhysicalB2State{T}) where {T}
  result = Dict{Tuple{CKPhysicalB2Kind,Int,Int},T}()
  for ((kind, harmonic, second_harmonic), value) in state.terms
    kind == CKPhysicalB2Bare || continue
    key = (CKPhysicalB2Solved, harmonic, second_harmonic)
    result[key] = get(result, key, state.zero_component) + value
  end
  return ck_physical_b2_state(result, state.zero_component)
end

function ck_physical_b2_terms(state::CKPhysicalB2State{T}) where {T}
  result = Dict{Tuple{Int,Int},T}()
  for ((kind, first_harmonic, second_harmonic), value) in state.terms
    kind == CKPhysicalB2PhysicalTwo || continue
    result[(first_harmonic, second_harmonic)] = value
  end
  return result
end

function ck_physical_b2_inverse_weight(mismatch::Int)
  iszero(mismatch) && return zero(CKPhysicalB2Exact)
  coefficient = ck_bernoulli_fourier_coefficient(ck_bernoulli_kernel(1), mismatch)
  return ck_physical_b2_im * get(coefficient, 0, zero(CKPhysicalB2Exact))
end

function ck_physical_b2_ordered_sideband_coefficient(
  state::CKPhysicalB2State{T}, first_sideband::Int, second_sideband::Int
) where {T}
  result = state.zero_component
  for ((first_harmonic, second_harmonic), value) in ck_physical_b2_terms(state)
    mismatch = first_harmonic - first_sideband
    iszero(mismatch) && continue
    second_harmonic + mismatch == second_sideband || continue
    result += ck_physical_b2_inverse_weight(mismatch) * value
  end
  return result
end

function ck_physical_b2_bosonic_sideband_coefficient(
  state::CKPhysicalB2State, first_sideband::Int, second_sideband::Int
)
  result = ck_physical_b2_ordered_sideband_coefficient(
    state, first_sideband, second_sideband
  )
  first_sideband == second_sideband && return result
  return result +
         ck_physical_b2_ordered_sideband_coefficient(state, second_sideband, first_sideband)
end

function ck_physical_b2_commutator_oracle(
  amplitudes::AbstractDict{Int,T},
  first_sideband::Int,
  second_sideband::Int,
  zero_component::T,
) where {T}
  result = zero_component
  sideband_assignments = if first_sideband == second_sideband
    ((first_sideband, second_sideband),)
  else
    ((first_sideband, second_sideband), (second_sideband, first_sideband))
  end

  for (first_harmonic, first_value) in amplitudes,
    (second_harmonic, second_value) in amplitudes,
    (left_sideband, right_sideband) in sideband_assignments

    mismatch = first_harmonic - left_sideband
    mismatch > 0 || continue
    second_harmonic + mismatch == right_sideband || continue
    commutator = first_value * second_value - second_value * first_value
    result += (-ck_physical_b2_im / mismatch) * commutator
  end
  return result
end

@testset "physical B2 homological kernel is the finite periodic Bernoulli sawtooth" begin
  kernel = ck_bernoulli_scale(ck_physical_b2_im, ck_bernoulli_kernel(1))
  @test ck_bernoulli_term_dictionary(kernel) ==
    Dict((1, 0) => CKPhysicalB2Exact(-1 // 2), (0, 1) => one(CKPhysicalB2Exact))
  @test isempty(ck_bernoulli_average(kernel))

  for mismatch in vcat(collect(-16:-1), collect(1:16))
    @test ck_physical_b2_inverse_weight(mismatch) == ck_physical_b2_im / mismatch
  end
end

@testset "order-only Bloch recurrence produces finite physical B2 kernel" begin
  zero_component = zeros(CKPhysicalB2Exact, 2, 2)
  identity_component = Matrix{CKPhysicalB2Exact}(I, 2, 2)
  amplitudes = Dict(
    -1 => CKPhysicalB2Exact[1 1; 0 -1],
    0 => CKPhysicalB2Exact[0 1; 1 1],
    1 => CKPhysicalB2Exact[2 -1; 1 0],
  )

  A1 = ck_physical_b2_bare(amplitudes, zero_component)
  identity_state = ck_physical_b2_identity(zero_component, identity_component)
  zero_state = CKPhysicalB2State(
    Dict{Tuple{CKPhysicalB2Kind,Int,Int},Matrix{CKPhysicalB2Exact}}(), zero_component
  )

  result = bloch_order_recurrence(
    [A1],
    2,
    identity_state,
    zero_state;
    product=ck_physical_b2_product,
    project_model=ck_physical_b2_project_model,
    solve_complement=ck_physical_b2_solve_complement,
  )

  @test result.products == 3
  @test count(key -> first(key) == CKPhysicalB2Model, keys(result.effective[1].terms)) == 3
  @test count(key -> first(key) == CKPhysicalB2Solved, keys(result.wave[1].terms)) == 3

  physical_terms = ck_physical_b2_terms(result.effective[2])
  @test length(physical_terms) == 9
  for (first_harmonic, first_value) in amplitudes,
    (second_harmonic, second_value) in amplitudes

    @test physical_terms[(first_harmonic, second_harmonic)] == second_value * first_value
  end

  @test all(key -> first(key) != CKPhysicalB2ComplementTwo, keys(result.effective[2].terms))

  for first_sideband in -8:8, second_sideband in first_sideband:8
    from_kernel = ck_physical_b2_bosonic_sideband_coefficient(
      result.effective[2], first_sideband, second_sideband
    )
    from_commutator = ck_physical_b2_commutator_oracle(
      amplitudes, first_sideband, second_sideband, zero_component
    )
    @test from_kernel == from_commutator
  end

  far_kernel = ck_physical_b2_bosonic_sideband_coefficient(result.effective[2], -32, 32)
  far_oracle = ck_physical_b2_commutator_oracle(amplitudes, -32, 32, zero_component)
  @test far_kernel == far_oracle
  @test !all(iszero, far_kernel)
end
