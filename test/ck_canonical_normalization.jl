using Test
using FloquetExpansions
using LinearAlgebra: I

const CKCanonicalExact = Complex{Rational{Int}}
const ck_canonical_im = CKCanonicalExact(0 // 1, 1 // 1)

function ck_canonical_inverse_weight(mismatch::Int)
  iszero(mismatch) && throw(ArgumentError("homological inverse requires nonzero mismatch"))
  return ck_canonical_im * (1 // mismatch)
end

function ck_canonical_endpoint_query(kernel, sidebands::Vector{Int})
  return FloquetExpansions.ck_endpoint_ordered_sideband_coefficient(
    kernel,
    fill(1, length(sidebands)),
    sidebands;
    inverse_weight=ck_canonical_inverse_weight,
  )
end

function ck_canonical_zero_sideband_query(kernel, outputs::Int)
  return ck_canonical_endpoint_query(kernel, zeros(Int, outputs))
end

function ck_canonical_grades(kernel)
  return sort!(
    unique(FloquetExpansions.ck_endpoint_output_number(key) for key in keys(kernel.terms))
  )
end

function ck_canonical_fixture()
  H0 = CKCanonicalExact[2 1 + ck_canonical_im; 1 - ck_canonical_im -1]
  H1 = CKCanonicalExact[1 2 - ck_canonical_im; -1 1 + ck_canonical_im]
  hamiltonian = Dict(0 => H0, 1 => H1, -1 => Matrix(adjoint(H1)))
  jumps = Dict(
    -1 => CKCanonicalExact[1 0; 2 ck_canonical_im],
    0 => CKCanonicalExact[0 1; -1 2],
    1 => CKCanonicalExact[1 - ck_canonical_im 2; 0 -1],
  )

  zero_component = zeros(CKCanonicalExact, 2, 2)
  identity_component = Matrix{CKCanonicalExact}(I, 2, 2)
  A1 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_jump_vertex(1, harmonic) => value for (harmonic, value) in jumps
    ),
    zero_component,
  )
  A2 = FloquetExpansions.ck_kernel_generator(
    Dict(
      FloquetExpansions.ck_drift_vertex(harmonic) => -ck_canonical_im * value for
      (harmonic, value) in hamiltonian
    ),
    zero_component,
  )
  identity_state = FloquetExpansions.ck_kernel_identity(
    0, identity_component, zero_component
  )
  zero_state = FloquetExpansions.ck_kernel_zero(0, zero_component)
  operations = FloquetExpansions.BlochOrderOperations(
    FloquetExpansions.ck_kernel_product,
    FloquetExpansions.ck_kernel_project_model,
    FloquetExpansions.ck_kernel_solve_complement,
  )
  return (; A1, A2, identity_state, zero_state, zero_component, operations)
end

struct CKCanonicalReference{T}
  components::Dict{Tuple{Int,Int},T}
  zero_component::T
end

function ck_reference_iszero(value)
  return all(iszero, value)
end

function ck_reference(components::Dict{Tuple{Int,Int},T}, zero_component::T) where {T}
  cleaned = Dict{Tuple{Int,Int},T}()
  for (key, value) in components
    ck_reference_iszero(value) || (cleaned[key] = value)
  end
  return CKCanonicalReference(cleaned, zero_component)
end

function Base.zero(value::CKCanonicalReference{T}) where {T}
  return CKCanonicalReference(Dict{Tuple{Int,Int},T}(), zero(value.zero_component))
end

function Base.getindex(value::CKCanonicalReference, harmonic::Int, grade::Int)
  return get(value.components, (harmonic, grade), value.zero_component)
end

function Base.:+(left::CKCanonicalReference{T}, right::CKCanonicalReference{T}) where {T}
  out = copy(left.components)
  for (key, value) in right.components
    out[key] = get(out, key, left.zero_component) + value
  end
  return ck_reference(out, left.zero_component)
end

function ck_reference_scale(value::CKCanonicalReference{T}, factor) where {T}
  return ck_reference(
    Dict(key => factor * component for (key, component) in value.components),
    value.zero_component,
  )
end

function ck_reference_product(
  left::CKCanonicalReference{T}, right::CKCanonicalReference{T}
) where {T}
  out = Dict{Tuple{Int,Int},T}()
  for ((left_harmonic, left_grade), left_component) in left.components
    for ((right_harmonic, right_grade), right_component) in right.components
      key = (left_harmonic + right_harmonic, left_grade + right_grade)
      out[key] = get(out, key, left.zero_component) + left_component * right_component
    end
  end
  return ck_reference(out, left.zero_component)
end

function ck_reference_commutator(left::CKCanonicalReference, right::CKCanonicalReference)
  return ck_reference_product(left, right) +
         ck_reference_scale(ck_reference_product(right, left), -1)
end

function ck_reference_project(value::CKCanonicalReference{T}) where {T}
  return ck_reference(
    Dict(key => component for (key, component) in value.components if iszero(first(key))),
    value.zero_component,
  )
end

function ck_reference_q_inverse(value::CKCanonicalReference{T}) where {T}
  out = Dict{Tuple{Int,Int},T}()
  for ((harmonic, grade), component) in value.components
    iszero(harmonic) && continue
    out[(harmonic, grade)] = ck_canonical_im * (1 // harmonic) * component
  end
  return ck_reference(out, value.zero_component)
end

function ck_reference_derivative(value::CKCanonicalReference{T}) where {T}
  out = Dict{Tuple{Int,Int},T}()
  for ((harmonic, grade), component) in value.components
    iszero(harmonic) && continue
    out[(harmonic, grade)] = (-ck_canonical_im * harmonic) * component
  end
  return ck_reference(out, value.zero_component)
end

function ck_reference_series_ad(
  generator::Vector{P}, values::Vector{P}
) where {P<:CKCanonicalReference}
  order = length(values)
  result = P[zero(first(values)) for _ in 1:order]
  for n in 1:order
    coefficient = zero(first(values))
    for j in 1:(n - 1)
      coefficient += ck_reference_commutator(generator[j], values[n - j])
    end
    result[n] = coefficient
  end
  return result
end

function ck_reference_hori_deprit(
  amplitude_orders::Vector{P}, order::Int
) where {P<:CKCanonicalReference}
  template = first(amplitude_orders)
  amplitude = P[zero(template) for _ in 1:order]
  for r in 1:min(order, length(amplitude_orders))
    amplitude[r] = amplitude_orders[r]
  end

  generator = P[zero(template) for _ in 1:order]
  effective = P[zero(template) for _ in 1:order]
  for n in 1:order
    source = zero(template)

    ad_amplitude = copy(amplitude)
    for power in 0:(n - 1)
      source += ck_reference_scale(ad_amplitude[n], (-1)^power // factorial(power))
      ad_amplitude = ck_reference_series_ad(generator, ad_amplitude)
    end

    derivative = P[ck_reference_derivative(component) for component in generator]
    ad_derivative = derivative
    for power in 0:(n - 1)
      weight = (-1)^(power + 1) // factorial(power + 1)
      source += ck_reference_scale(ad_derivative[n], weight)
      ad_derivative = ck_reference_series_ad(generator, ad_derivative)
    end

    effective[n] = ck_reference_project(source)
    generator[n] = ck_reference_q_inverse(source)
  end
  return effective
end

function ck_canonical_reference_fixture()
  H0 = CKCanonicalExact[2 1 + ck_canonical_im; 1 - ck_canonical_im -1]
  H1 = CKCanonicalExact[1 2 - ck_canonical_im; -1 1 + ck_canonical_im]
  hamiltonian = Dict(0 => H0, 1 => H1, -1 => Matrix(adjoint(H1)))
  jumps = Dict(
    -1 => CKCanonicalExact[1 0; 2 ck_canonical_im],
    0 => CKCanonicalExact[0 1; -1 2],
    1 => CKCanonicalExact[1 - ck_canonical_im 2; 0 -1],
  )
  zero_component = zeros(CKCanonicalExact, 2, 2)
  A1 = ck_reference(
    Dict((harmonic, 1) => value for (harmonic, value) in jumps), zero_component
  )
  A2 = ck_reference(
    Dict((harmonic, 0) => -ck_canonical_im * value for (harmonic, value) in hamiltonian),
    zero_component,
  )
  return ck_reference_hori_deprit([A1, A2], 5)
end

@testset "CK endpoint canonical algebra preserves block-local complement constraints" begin
  qplus_physical = FloquetExpansions.ck_kernel_solve_complement(
    FloquetExpansions.ck_kernel_generator(
      Dict(FloquetExpansions.ck_drift_vertex(1) => CKCanonicalExact(1)), CKCanonicalExact(0)
    ),
  )
  qminus_physical = FloquetExpansions.ck_kernel_solve_complement(
    FloquetExpansions.ck_kernel_generator(
      Dict(FloquetExpansions.ck_drift_vertex(-1) => CKCanonicalExact(1)),
      CKCanonicalExact(0),
    ),
  )
  qplus = FloquetExpansions.ck_endpoint_kernel(qplus_physical)
  qminus = FloquetExpansions.ck_endpoint_kernel(qminus_physical)
  product = FloquetExpansions.ck_endpoint_product(qplus, qminus)
  projected = FloquetExpansions.ck_endpoint_project_model(product)

  @test isempty(FloquetExpansions.ck_kernel_product(qplus_physical, qminus_physical).terms)
  @test length(projected.terms) == 1
  key = first(keys(projected.terms))
  @test key.vertices ==
    [FloquetExpansions.ck_drift_vertex(-1), FloquetExpansions.ck_drift_vertex(1)]
  @test key.resolvent_intervals == [
    FloquetExpansions.CKEndpointConstraint(0, 1),
    FloquetExpansions.CKEndpointConstraint(1, 2),
  ]
  @test key.model_intervals == [FloquetExpansions.CKEndpointConstraint(0, 2)]
  @test ck_canonical_zero_sideband_query(projected, 0) == CKCanonicalExact(1)

  identity = FloquetExpansions.ck_endpoint_kernel(
    FloquetExpansions.ck_kernel_identity(0, CKCanonicalExact(1), CKCanonicalExact(0))
  )
  @test FloquetExpansions.ck_endpoint_product(identity, qplus) == qplus
  @test FloquetExpansions.ck_endpoint_product(qplus, identity) == qplus
end

@testset "CK canonical normalization reproduces lower-order gauge identities" begin
  fixture = ck_canonical_fixture()
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    5,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective, recurrence.wave, 5, fixture.identity_state, fixture.zero_state
  )

  @test isempty(canonical.static_factor[1].terms)
  X1 = FloquetExpansions.ck_endpoint_kernel(recurrence.wave[1])
  expected_N2 = FloquetExpansions.ck_endpoint_scale(
    FloquetExpansions.ck_endpoint_project_model(
      FloquetExpansions.ck_endpoint_product(X1, X1)
    ),
    1 // 2,
  )
  @test canonical.static_factor[2] == expected_N2

  B1 = FloquetExpansions.ck_endpoint_kernel(recurrence.effective[1])
  B3 = FloquetExpansions.ck_endpoint_kernel(recurrence.effective[3])
  commutator_B1_N2 =
    FloquetExpansions.ck_endpoint_product(B1, canonical.static_factor[2]) -
    FloquetExpansions.ck_endpoint_product(canonical.static_factor[2], B1)
  @test canonical.effective[3] == B3 + commutator_B1_N2

  for sideband in -4:4
    @test ck_canonical_endpoint_query(canonical.effective[3], [sideband]) ==
      FloquetExpansions.ck_kernel_ordered_sideband_coefficient(
      recurrence.effective[3], [1], [sideband]; inverse_weight=ck_canonical_inverse_weight
    )
  end

  zero_endpoint = FloquetExpansions.ck_endpoint_kernel(fixture.zero_state)
  for coefficient in canonical.log_embedding
    @test FloquetExpansions.ck_endpoint_project_model(coefficient) == zero_endpoint
  end
end

@testset "CK canonical normalization matches independent Hori-Deprit through order five" begin
  fixture = ck_canonical_fixture()
  recurrence = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    5,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence.effective, recurrence.wave, 5, fixture.identity_state, fixture.zero_state
  )
  reference = ck_canonical_reference_fixture()

  for order in 1:5, grade in 0:order
    @test ck_canonical_zero_sideband_query(canonical.effective[order], grade) ==
      reference[order][0, grade]
  end

  @test ck_canonical_grades(canonical.static_factor[2]) == [2]
  @test ck_canonical_grades(canonical.static_factor[3]) == [1, 3]
  @test ck_canonical_grades(canonical.static_factor[4]) == [0, 2, 4]
  @test !iszero(ck_canonical_zero_sideband_query(canonical.static_factor[3], 1))
  @test !iszero(ck_canonical_zero_sideband_query(canonical.static_factor[4], 0))

  raw_B5 = FloquetExpansions.ck_endpoint_kernel(recurrence.effective[5])
  raw_grade1 = ck_canonical_zero_sideband_query(raw_B5, 1)
  canonical_grade1 = ck_canonical_zero_sideband_query(canonical.effective[5], 1)
  @test raw_grade1 != canonical_grade1

  B1 = FloquetExpansions.ck_endpoint_kernel(recurrence.effective[1])
  B2 = FloquetExpansions.ck_endpoint_kernel(recurrence.effective[2])
  N3 = canonical.static_factor[3]
  N4 = canonical.static_factor[4]
  correction =
    FloquetExpansions.ck_endpoint_product(B1, N4) -
    FloquetExpansions.ck_endpoint_product(N4, B1) +
    FloquetExpansions.ck_endpoint_product(B2, N3) -
    FloquetExpansions.ck_endpoint_product(N3, B2)
  @test raw_grade1 + ck_canonical_zero_sideband_query(correction, 1) == canonical_grade1
end

@testset "CK canonical normalization preserves prefixes and period reconstruction" begin
  fixture = ck_canonical_fixture()
  recurrence4 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    4,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  recurrence5 = FloquetExpansions.evaluate_bloch_order_recurrence(
    [fixture.A1, fixture.A2],
    5,
    fixture.identity_state,
    fixture.zero_state,
    fixture.operations,
  )
  canonical4 = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence4.effective, recurrence4.wave, 4, fixture.identity_state, fixture.zero_state
  )
  canonical5 = FloquetExpansions.evaluate_ck_canonical_normalization(
    recurrence5.effective, recurrence5.wave, 5, fixture.identity_state, fixture.zero_state
  )

  @test canonical5.static_factor[1:3] == canonical4.static_factor
  @test canonical5.inverse_static_factor[1:3] == canonical4.inverse_static_factor
  @test canonical5.wave[1:3] == canonical4.wave
  @test canonical5.effective[1:4] == canonical4.effective

  raw_period = FloquetExpansions.evaluate_ck_period_amplitude(
    recurrence4.effective, recurrence4.wave, 4, fixture.identity_state, fixture.zero_state
  )
  endpoint_period = FloquetExpansions.evaluate_ck_period_amplitude(
    [FloquetExpansions.ck_endpoint_kernel(value) for value in recurrence4.effective],
    [FloquetExpansions.ck_endpoint_kernel(value) for value in recurrence4.wave],
    4,
    FloquetExpansions.ck_endpoint_kernel(fixture.identity_state),
    FloquetExpansions.ck_endpoint_kernel(fixture.zero_state),
  )
  @test endpoint_period.amplitude == raw_period.amplitude
  @test endpoint_period.endpoint_wave == raw_period.endpoint_wave
  @test endpoint_period.endpoint_wave_inverse == raw_period.endpoint_wave_inverse
  @test endpoint_period.slow_propagator == raw_period.slow_propagator

  @test FloquetExpansions.ck_canonical_total_products(canonical5.counts) > 0
end
