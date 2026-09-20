using Test
using LinearAlgebra: I

const GaugeExact = Complex{Rational{Int}}
const gauge_im = GaugeExact(0 // 1, 1 // 1)

struct GaugeOpenLeg{T}
  components::Dict{Tuple{Int,Int},T}
  zero_component::T
end

gauge_component_iszero(value) = iszero(value)
gauge_component_iszero(value::AbstractMatrix) = all(iszero, value)

function gauge_open_leg(components::Dict{Tuple{Int,Int},T}, zero_component::T) where {T}
  cleaned = Dict{Tuple{Int,Int},T}()
  for (key, value) in components
    gauge_component_iszero(value) || (cleaned[key] = value)
  end
  return GaugeOpenLeg(cleaned, zero_component)
end

function Base.zero(value::GaugeOpenLeg{T}) where {T}
  return GaugeOpenLeg(Dict{Tuple{Int,Int},T}(), zero(value.zero_component))
end

function Base.getindex(value::GaugeOpenLeg, harmonic::Int, grade::Int)
  return get(value.components, (harmonic, grade), value.zero_component)
end

function Base.:(==)(left::GaugeOpenLeg, right::GaugeOpenLeg)
  all_keys = union(keys(left.components), keys(right.components))
  return all(left[key...] == right[key...] for key in all_keys)
end

function Base.:+(left::GaugeOpenLeg{T}, right::GaugeOpenLeg{T}) where {T}
  out = copy(left.components)
  for (key, value) in right.components
    out[key] = get(out, key, left.zero_component) + value
  end
  return gauge_open_leg(out, left.zero_component)
end

function Base.:-(left::GaugeOpenLeg{T}, right::GaugeOpenLeg{T}) where {T}
  out = copy(left.components)
  for (key, value) in right.components
    out[key] = get(out, key, left.zero_component) - value
  end
  return gauge_open_leg(out, left.zero_component)
end

function Base.:-(value::GaugeOpenLeg{T}) where {T}
  return gauge_open_leg(
    Dict(key => -component for (key, component) in value.components), value.zero_component
  )
end

function Base.:*(scalar::Number, value::GaugeOpenLeg{T}) where {T}
  return gauge_open_leg(
    Dict(key => scalar * component for (key, component) in value.components),
    value.zero_component,
  )
end

Base.:*(value::GaugeOpenLeg, scalar::Number) = scalar * value

function gauge_product(left::GaugeOpenLeg{T}, right::GaugeOpenLeg{T}) where {T}
  out = Dict{Tuple{Int,Int},T}()
  for ((left_harmonic, left_grade), left_component) in left.components
    for ((right_harmonic, right_grade), right_component) in right.components
      key = (left_harmonic + right_harmonic, left_grade + right_grade)
      term = left_component * right_component
      out[key] = get(out, key, left.zero_component) + term
    end
  end
  return gauge_open_leg(out, left.zero_component)
end

function gauge_commutator(left::GaugeOpenLeg, right::GaugeOpenLeg)
  return gauge_product(left, right) - gauge_product(right, left)
end

function gauge_project(value::GaugeOpenLeg{T}) where {T}
  return gauge_open_leg(
    Dict(key => component for (key, component) in value.components if iszero(first(key))),
    value.zero_component,
  )
end

function gauge_q_inverse(value::GaugeOpenLeg{T}) where {T}
  out = Dict{Tuple{Int,Int},T}()
  for ((harmonic, grade), component) in value.components
    iszero(harmonic) && continue
    out[(harmonic, grade)] = gauge_im * (1 // harmonic) * component
  end
  return gauge_open_leg(out, value.zero_component)
end

function gauge_derivative(value::GaugeOpenLeg{T}) where {T}
  out = Dict{Tuple{Int,Int},T}()
  for ((harmonic, grade), component) in value.components
    iszero(harmonic) && continue
    out[(harmonic, grade)] = (-gauge_im * harmonic) * component
  end
  return gauge_open_leg(out, value.zero_component)
end

function gauge_grade(value::GaugeOpenLeg{T}, grade::Int) where {T}
  return gauge_open_leg(
    Dict(key => component for (key, component) in value.components if last(key) == grade),
    value.zero_component,
  )
end

function gauge_grades(value::GaugeOpenLeg)
  return sort!(unique(last(key) for key in keys(value.components)))
end

function gauge_identity(template::GaugeOpenLeg{T}, identity_component::T) where {T}
  return gauge_open_leg(Dict((0, 0) => identity_component), template.zero_component)
end

struct GaugeBlochReference{P}
  wave::Vector{P}
  effective::Vector{P}
end

function gauge_bloch_reference(
  amplitude_orders::Vector{P}, order::Int, identity_component
) where {P<:GaugeOpenLeg}
  template = first(amplitude_orders)
  identity = gauge_identity(template, identity_component)
  wave = P[]
  effective = P[]

  for n in 1:order
    residual = zero(template)
    for r in 1:min(n, length(amplitude_orders))
      previous = n == r ? identity : wave[n - r]
      residual = residual + gauge_product(amplitude_orders[r], previous)
    end
    for j in 1:(n - 1)
      residual = residual - gauge_product(wave[j], effective[n - j])
    end
    push!(effective, gauge_project(residual))
    push!(wave, gauge_q_inverse(residual))
  end

  return GaugeBlochReference(wave, effective)
end

struct GaugeHoriDepritReference{P}
  generator::Vector{P}
  effective::Vector{P}
end

function gauge_series_ad(generator::Vector{P}, values::Vector{P}) where {P<:GaugeOpenLeg}
  order = length(values)
  result = [zero(first(values)) for _ in 1:order]
  for n in 1:order
    coefficient = zero(first(values))
    for j in 1:(n - 1)
      coefficient = coefficient + gauge_commutator(generator[j], values[n - j])
    end
    result[n] = coefficient
  end
  return result
end

function gauge_hori_deprit_reference(
  amplitude_orders::Vector{P}, order::Int
) where {P<:GaugeOpenLeg}
  template = first(amplitude_orders)
  amplitude = [zero(template) for _ in 1:order]
  for r in 1:min(order, length(amplitude_orders))
    amplitude[r] = amplitude_orders[r]
  end

  generator = [zero(template) for _ in 1:order]
  effective = [zero(template) for _ in 1:order]

  for n in 1:order
    source = zero(template)

    ad_amplitude = copy(amplitude)
    for power in 0:(n - 1)
      source = source + ((-1)^power // factorial(power)) * ad_amplitude[n]
      ad_amplitude = gauge_series_ad(generator, ad_amplitude)
    end

    derivative = [gauge_derivative(component) for component in generator]
    ad_derivative = derivative
    for power in 0:(n - 1)
      weight = (-1)^(power + 1) // factorial(power + 1)
      source = source + weight * ad_derivative[n]
      ad_derivative = gauge_series_ad(generator, ad_derivative)
    end

    effective[n] = gauge_project(source)
    generator[n] = gauge_q_inverse(source)
  end

  return GaugeHoriDepritReference(generator, effective)
end

function gauge_positive_series_product(
  left::Vector{P}, right::Vector{P}
) where {P<:GaugeOpenLeg}
  order = length(left)
  result = [zero(first(left)) for _ in 1:order]
  for n in 1:order
    coefficient = zero(first(left))
    for j in 1:(n - 1)
      coefficient = coefficient + gauge_product(left[j], right[n - j])
    end
    result[n] = coefficient
  end
  return result
end

function gauge_series_log(values::Vector{P}) where {P<:GaugeOpenLeg}
  order = length(values)
  result = [zero(first(values)) for _ in 1:order]
  powers = copy(values)
  for power in 1:order
    weight = (-1)^(power + 1) // power
    for n in 1:order
      result[n] = result[n] + weight * powers[n]
    end
    powers = gauge_positive_series_product(values, powers)
  end
  return result
end

function gauge_static_factor(wave::Vector{P}, order::Int) where {P<:GaugeOpenLeg}
  order >= 1 || throw(ArgumentError("order must be positive"))
  length(wave) >= order || throw(ArgumentError("insufficient Bloch wave data"))

  template = first(wave)
  static = [zero(template) for _ in 1:order]
  normalized = [zero(template) for _ in 1:order]

  for n in 1:order
    prefactor = wave[n]
    for j in 1:(n - 1)
      prefactor = prefactor + gauge_product(wave[j], static[n - j])
    end

    candidate = copy(normalized)
    candidate[n] = prefactor
    candidate_log = gauge_series_log(candidate)
    static[n] = -gauge_project(candidate_log[n])
    normalized[n] = prefactor + static[n]
  end

  return (; static, normalized, log=gauge_series_log(normalized))
end

function gauge_static_inverse(static::Vector{P}) where {P<:GaugeOpenLeg}
  order = length(static)
  inverse = [zero(first(static)) for _ in 1:order]
  for n in 1:order
    coefficient = static[n]
    for j in 1:(n - 1)
      coefficient = coefficient + gauge_product(static[j], inverse[n - j])
    end
    inverse[n] = -coefficient
  end
  return inverse
end

function gauge_canonical_effective(
  effective::Vector{P}, static::Vector{P}, order::Int
) where {P<:GaugeOpenLeg}
  inverse = gauge_static_inverse(static)
  right = [zero(first(effective)) for _ in 1:order]
  canonical = [zero(first(effective)) for _ in 1:order]

  for n in 1:order
    coefficient = effective[n]
    for j in 1:(n - 1)
      coefficient = coefficient + gauge_product(effective[j], static[n - j])
    end
    right[n] = coefficient
  end

  for n in 1:order
    coefficient = right[n]
    for j in 1:(n - 1)
      coefficient = coefficient + gauge_product(inverse[j], right[n - j])
    end
    canonical[n] = coefficient
  end
  return canonical
end

function gauge_fixture()
  H0 = GaugeExact[2 1 + gauge_im; 1 - gauge_im -1]
  H1 = GaugeExact[1 2 - gauge_im; -1 1 + gauge_im]
  hamiltonian = Dict(0 => H0, 1 => H1, -1 => Matrix(adjoint(H1)))

  jumps = Dict(
    -1 => GaugeExact[1 0; 2 gauge_im],
    0 => GaugeExact[0 1; -1 2],
    1 => GaugeExact[1 - gauge_im 2; 0 -1],
  )

  zero_component = zeros(GaugeExact, 2, 2)
  identity_component = Matrix{GaugeExact}(I, 2, 2)
  A1 = gauge_open_leg(
    Dict((harmonic, 1) => value for (harmonic, value) in jumps), zero_component
  )
  A2 = gauge_open_leg(
    Dict((harmonic, 0) => -gauge_im * value for (harmonic, value) in hamiltonian),
    zero_component,
  )
  return (; A1, A2, identity_component, zero_component)
end

@testset "open-leg generic Hori-Deprit recurrence reproduces the cubic oracle" begin
  fixture = gauge_fixture()
  hd = gauge_hori_deprit_reference([fixture.A1, fixture.A2], 5)

  G1 = gauge_q_inverse(fixture.A1)
  derivative_G1 = gauge_derivative(G1)
  F2 =
    fixture.A2 - gauge_commutator(G1, fixture.A1) +
    (1 // 2) * gauge_commutator(G1, derivative_G1)
  B2 = gauge_project(F2)
  G2 = gauge_q_inverse(F2)

  F3 =
    -gauge_commutator(G1, fixture.A2) - gauge_commutator(G2, fixture.A1) +
    (1 // 2) * gauge_commutator(G1, gauge_commutator(G1, fixture.A1)) +
    (1 // 2) * gauge_commutator(G1, gauge_derivative(G2)) +
    (1 // 2) * gauge_commutator(G2, derivative_G1) -
    (1 // 6) * gauge_commutator(G1, gauge_commutator(G1, derivative_G1))

  @test hd.generator[1] == G1
  @test hd.effective[2] == B2
  @test hd.generator[2] == G2
  @test hd.effective[3] == gauge_project(F3)
end

@testset "open-leg Bloch canonical normalization matches Hori-Deprit through order five" begin
  fixture = gauge_fixture()
  amplitude_orders = [fixture.A1, fixture.A2]
  bloch = gauge_bloch_reference(amplitude_orders, 5, fixture.identity_component)
  hd = gauge_hori_deprit_reference(amplitude_orders, 5)

  @test bloch.effective[2] == hd.effective[2]
  @test gauge_grade(bloch.effective[3], 1) == gauge_grade(hd.effective[3], 1)
  @test bloch.effective[3] != hd.effective[3]

  normalization = gauge_static_factor(bloch.wave, 4)
  static = normalization.static
  canonical = gauge_canonical_effective(bloch.effective, static, 5)

  @test all(gauge_project(component) == zero(component) for component in normalization.log)
  @test static[1] == zero(static[1])
  @test static[2] == (1 // 2) * gauge_project(gauge_product(bloch.wave[1], bloch.wave[1]))
  @test gauge_grades(static[2]) == [2]
  @test gauge_grades(static[3]) == [1, 3]
  @test gauge_grades(static[4]) == [0, 2, 4]
  @test !gauge_component_iszero(static[4][0, 0])

  for order in 1:5
    @test canonical[order] == hd.effective[order]
  end

  raw_grade1 = gauge_grade(bloch.effective[5], 1)
  canonical_grade1 = gauge_grade(canonical[5], 1)
  hd_grade1 = gauge_grade(hd.effective[5], 1)
  @test raw_grade1 != hd_grade1

  grade1_correction =
    gauge_commutator(gauge_grade(bloch.effective[1], 1), gauge_grade(static[4], 0)) +
    gauge_commutator(gauge_grade(bloch.effective[2], 0), gauge_grade(static[3], 1))
  @test raw_grade1 + grade1_correction == canonical_grade1
  @test canonical_grade1 == hd_grade1
end

@testset "open-leg order-five references preserve lower-order prefixes" begin
  fixture = gauge_fixture()
  amplitude_orders = [fixture.A1, fixture.A2]

  bloch4 = gauge_bloch_reference(amplitude_orders, 4, fixture.identity_component)
  bloch5 = gauge_bloch_reference(amplitude_orders, 5, fixture.identity_component)
  hd4 = gauge_hori_deprit_reference(amplitude_orders, 4)
  hd5 = gauge_hori_deprit_reference(amplitude_orders, 5)
  normalization4 = gauge_static_factor(bloch4.wave, 3)
  normalization5 = gauge_static_factor(bloch5.wave, 4)

  @test bloch5.effective[1:4] == bloch4.effective
  @test bloch5.wave[1:4] == bloch4.wave
  @test hd5.effective[1:4] == hd4.effective
  @test hd5.generator[1:4] == hd4.generator
  @test normalization5.static[1:3] == normalization4.static
end