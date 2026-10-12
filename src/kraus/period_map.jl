"""
    JumpAmplitude

One continuous Kraus operator of a [`PeriodMap`](@ref): the retained amplitude ``K(u)`` of the
jump records with ``q`` jumps at the fractional times ``0 < u_1 < \\dots < u_q < 1`` of the period,
with its weight ``w``, nonnegative when the jump rates are. It contributes
``w \\int_{0<u_1<\\dots<u_q<1} K(u) \\rho K(u)^\\dagger \\, du`` to the map. Read it with
[`weight`](@ref), [`amplitude`](@ref) and [`jump_times`](@ref).
"""
struct JumpAmplitude
  weight::SQA.CNum
  amplitude::SQA.QAdd
  jump_times::Vector{Symbolics.Num}
end

"""
    weight(K::JumpAmplitude) -> CNum

Return the weight ``T^q \\gamma_{a_1} \\cdots \\gamma_{a_q}`` of the jump records of `K`,
with ``T = 2\\pi/\\omega``.
"""
weight(K::JumpAmplitude) = K.weight

"""
    amplitude(K::JumpAmplitude) -> QAdd

Return the operator ``K(u)``, a polynomial in ``1/\\omega`` whose coefficients are polynomials in
the [`jump_times`](@ref) ``u_j`` times ``e^{-2\\pi i k u_j}``.
"""
amplitude(K::JumpAmplitude) = K.amplitude

"""
    jump_times(K::JumpAmplitude) -> Vector{Num}

Return the fractional jump times ``u_1, \\dots, u_q`` of `K`, ordered as the jumps occur. The
vector is empty for the no-jump amplitude.
"""
jump_times(K::JumpAmplitude) = K.jump_times

"""
    PeriodMap

Result of [`period_map`](@ref): a completely positive map ``\\Phi`` that agrees with the exact
evolution from ``t = 0`` to ``t = T = 2\\pi/\\omega`` through ``\\omega^{-N}``, where ``N`` is
[`order`](@ref). Read it with [`superoperator`](@ref), [`jump_amplitudes`](@ref) and
[`trace_defect`](@ref) rather than by field access.
"""
struct PeriodMap
  superoperator::Liouvillian
  amplitudes::Vector{JumpAmplitude}
  trace_defect::SQA.QAdd
  order::Int
end

function Base.show(io::IO, ::MIME"text/plain", Φ::PeriodMap)
  return print(io, "PeriodMap of order ", Φ.order)
end
Base.show(io::IO, Φ::PeriodMap) = show(io, MIME"text/plain"(), Φ)

"""
    order(Φ::PeriodMap) -> Int

Return the inverse-frequency order ``N`` through which `Φ` agrees with the exact period map.
"""
order(Φ::PeriodMap) = Φ.order

"""
    superoperator(Φ::PeriodMap) -> Liouvillian

Return the period map as the linear map ``\\rho \\mapsto \\sum_j c_j A_j \\rho B_j``. Its
coefficients are exact polynomials in ``\\pi``, the parameters, and ``1/\\omega``. The map is
completely positive for every drive frequency and nonnegative jump rates, not only through the
retained order, because it is a sum over jump records of squared jump amplitudes; see
[Floquet–Kraus period map](@ref floquet-kraus-theory).

`SQA.simplify` evaluates ``\\pi`` in floating point, so simplify a coefficient only after
substituting values, or substitute with `fold=Val(true)`.
"""
superoperator(Φ::PeriodMap) = Φ.superoperator

"""
    trace_defect(Φ::PeriodMap) -> QAdd

Return ``\\Phi^\\dagger(1) - 1``, which measures the failure of trace preservation. It starts at
``\\omega^{-N-1}``, where ``N`` is [`order`](@ref); the construction checks exactly that every
lower order vanishes. The map ``\\rho \\mapsto \\Phi(S\\rho S)`` with ``S = (1 + D)^{-1/2}``, ``D``
the trace defect, is exactly trace preserving and still agrees with the exact period map through
``\\omega^{-N}``, wherever ``1 + D`` is positive definite.
"""
trace_defect(Φ::PeriodMap) = Φ.trace_defect

"""
    jump_amplitudes(Φ::PeriodMap) -> Vector{JumpAmplitude}

Return the period map in continuous Kraus form,
``\\Phi(\\rho) = \\sum_K w_K \\int K(u) \\rho K(u)^\\dagger \\, du``, with one
[`JumpAmplitude`](@ref) per sequence of jump channels. The integral over the jump times reproduces
[`superoperator`](@ref) exactly. A finite Kraus form follows at substituted values from a
factorization of the Gram matrix of the monomials in ``u``.
"""
jump_amplitudes(Φ::PeriodMap) = Φ.amplitudes

struct LoweredChannel{X,F}
  operator::Harmonics{X}
  rate::F
end

function lower_harmonics(lowering::SQALowering, P::PeriodicGenerator{SQA.QAdd})
  return Dict(k => lower_qadd(lowering, P[k]) for k in keys(P))
end

function check_time_independent_rate(rate::SQA.CNum, wd::Symbolics.Num, t::Symbolics.Num)
  found = Symbolics.get_variables(coefficient_number(rate))
  any(v -> isequal(v, Symbolics.unwrap(t)) || isequal(v, Symbolics.unwrap(wd)), found) &&
    throw(
      ArgumentError(
        "the period map needs time-independent jump rates, found $rate; put the time " *
        "dependence into the jump operator",
      ),
    )
  return rate
end

function lowered_channels(
  lowering::SQALowering{T}, provenance::MicroscopicProvenance, wd, t
) where {T}
  operators = [provenance.collapse_operators; provenance.jump_operators]
  rates = [
    fill(one(T), length(provenance.collapse_operators))
    [
      lower_coefficient(lowering, check_time_independent_rate(rate, wd, t)) for
      rate in provenance.jump_rates
    ]
  ]
  return LoweredChannel{AlgebraOperator{T,real(T)},T}[
    LoweredChannel(lower_harmonics(lowering, harmonics(operator, wd, t)), rate) for
    (operator, rate) in zip(operators, rates)
  ]
end

function no_jump_generator(
  H::Harmonics{X}, channels::Vector{<:LoweredChannel{X}}, zero_element::X
) where {X}
  Q = Harmonics{X}(k => -im * value for (k, value) in H)
  for channel in channels
    adjoint_operator = Harmonics{X}(-k => adjoint(value) for (k, value) in channel.operator)
    decay = harmonic_product(adjoint_operator, channel.operator)
    Q = harmonic_combine(Q, decay, 1, -channel.rate / 2)
  end
  for (k, value) in collect(Q)
    iszero(value) && delete!(Q, k)
  end
  haskey(Q, 0) || (Q[0] = zero_element)
  return Q
end

function period_map_lowering(generator::PeriodicGenerator{Liouvillian}, order::Int)
  parameters = checked_native_parameters(
    generator, MAXIMUM_PARAMETERS - 1, "the period map"
  )
  times = [Symbolics.variable(:u, j; T=Real) for j in 1:order]
  any(u -> any(p -> isequal(p, u), parameters), times) && throw(
    ArgumentError(
      "the jump times are named u₁, …, u_$order; rename the parameters that use these names",
    ),
  )
  T = Complex{ParameterScalar}
  values = T[parameter_variable(i) for i in eachindex(parameters)]
  lowering = SQALowering{T}(generator, parameters, values)
  coefficients = ParameterCoefficients(
    [parameters; Symbolics.Num(π)], ones(Int, length(parameters) + 1)
  )
  twopi = T(2 * parameter_variable(length(parameters) + 1))
  return lowering, coefficients, twopi, times
end

struct JumpRecords{X,F,L<:CoefficientLift}
  jumps::Vector{Vector{Harmonics{X}}}
  powers::Vector{Vector{X}}
  left::Vector{X}
  right::Vector{X}
  rates::Vector{F}
  twopi::F
  integrals::SimplexIntegrals{F}
  lift::L
  drive_scales::Vector{SQA.CNum}
  times::Vector{Symbolics.Num}
  time_values::Dict{TimeMonomial,SQA.CNum}
end

function record_rate(word::Vector{Int}, records::JumpRecords)
  rate = records.twopi^length(word)
  for channel in word
    rate *= records.rates[channel]
  end
  return rate
end

function record_map!(
  series::Vector, q::Int, kernel::TimeKernel{X}, rate, integrals::SimplexIntegrals
) where {X}
  adjoints = Dict(key => adjoint(B) for (key, B) in kernel.terms)
  for ((r, μ), A) in kernel.terms
    weighted = Dict{Int,X}()
    for ((s, ν), B) in adjoints
      weight = gram_entry(integrals, μ, ν)
      iszero(weight) || accumulate!(weighted, s, (rate * weight) * B)
    end
    for (s, B) in weighted
      series[q + r + s + 1] += sandwich(A, B)
    end
  end
  return series
end

function time_value(records::JumpRecords, monomial::TimeMonomial)
  return get!(() -> time_monomial(monomial, records.times), records.time_values, monomial)
end

function lift_amplitude(records::JumpRecords, kernel::TimeKernel, rate, q::Int)
  weight = lift_coefficient(records.lift.coefficient, rate) * records.drive_scales[q + 1]
  pieces = (
    (X, records.drive_scales[r + 1] * time_value(records, monomial)) for
    ((r, monomial), X) in kernel.terms
  )
  return JumpAmplitude(weight, lift_scaled_sum(records.lift, pieces), records.times[1:q])
end

function record!(series, amplitudes, records::JumpRecords, kernel, word, R::Int)
  amplitude = conjugate_kernel(kernel, records.left, records.right, R)
  isempty(amplitude.terms) && return series
  rate = record_rate(word, records)
  record_map!(series, length(word), amplitude, rate, records.integrals)
  push!(amplitudes, lift_amplitude(records, amplitude, rate, length(word)))
  return series
end

function jump_records!(
  series, amplitudes, records::JumpRecords, kernel, word::Vector{Int}, q::Int, R::Int
)
  isempty(kernel.terms) && return series
  j = length(word)
  j == q && return record!(series, amplitudes, records, kernel, word, R)
  for (channel, jump) in pairs(records.jumps)
    next = insert_jump(kernel, jump, j + 1, R)
    push!(word, channel)
    jump_records!(
      series,
      amplitudes,
      records,
      drift(next, records.powers, records.twopi, j + 2, j + 1, R),
      word,
      q,
      R,
    )
    pop!(word)
  end
  return series
end

function dual_identity(S::AlgebraSuperoperator{T,R}) where {T,R}
  terms = Dict{Monomial,T}()
  for ((left, right), c) in S.terms,
    (monomial, d) in monomial_product(S.algebra, right, left)

    accumulate!(terms, monomial, c * d)
  end
  return AlgebraOperator{T,R}(S.algebra, terms)
end

function graded_trace_defect(series::Vector{<:AlgebraSuperoperator}, identity, N::Int)
  defect = [dual_identity(S) for S in series]
  defect[1] -= identity
  all(iszero, defect[1:(N + 1)]) ||
    error("the period map failed its trace-preservation certificate")
  return defect
end

function time_monomial((phases, powers)::TimeMonomial, times::Vector{Symbolics.Num})
  value = Complex{Symbolics.Num}(1)
  for (u, k, p) in zip(times, phases, powers)
    angle = 2 * Symbolics.Num(π) * k * u
    value *= u^p * Complex(cos(angle), -sin(angle))
  end
  return convert(SQA.CNum, value)
end

function lift_graded_superoperator(lift::CoefficientLift, series::Vector, scales)
  result = Liouvillian(LiouvillianTerms())
  for (S, scale) in zip(series, scales)
    iszero(S) || lift_superoperator!(result, lift, S, scale)
  end
  return result
end

"""
    period_map(H::QField, ωd, t, order; channels=(), algorithm=HoriDeprit()) -> PeriodMap

Expand the one-period map of the periodically driven Lindblad system

```math
\\dot\\rho = -i[H(t), \\rho] + \\sum_a \\gamma_a \\mathcal{D}[L_a(t)]\\rho
```

from ``t = 0`` to ``t = T = 2\\pi/\\omega_d`` in powers of ``1/\\omega_d``, as a completely
positive map. The result agrees with the exact period map through
``\\omega_d^{-\\text{order}}``, so its logarithm divided by ``T`` agrees with the stroboscopic
effective generator through ``\\omega_d^{1-\\text{order}}``, as for
[`floquet_expansion`](@ref) with the same `order`.

The map is never assumed to have a time-independent GKSL logarithm. It is the sum over
quantum-jump records of the period of squared jump amplitudes, integrated exactly over the jump
times. Every amplitude is retained through the requested order and squared without truncation, so
the map is completely positive for every drive frequency, and its [`trace_defect`](@ref) starts
beyond the retained order. A symbolic rate enters as the nonnegative quantity that [`jump`](@ref)
assumes it to be, and complete positivity holds only where that assumption does.

`algorithm` selects the [`ExpansionAlgorithm`](@ref) for the Floquet frame of the no-jump
evolution: [`HoriDeprit`](@ref), the default, or [`BlochFeshbach`](@ref). Both give the same map.
The construction is described in [Floquet–Kraus period map](@ref floquet-kraus-theory).

# Arguments

- `H`: Symbolic time-dependent Hamiltonian, Hermitian at every time.
- `ωd`: Symbolic drive frequency.
- `t`: Symbolic time variable.
- `order`: Inverse-frequency order ``N \\ge 1`` through which the map is exact.
- `channels`: Tuple or vector of [`collapse`](@ref) and [`jump`](@ref) values.

# Restrictions

- Jump rates must be time independent; a periodic rate belongs in the operator.
- Coefficients must be exact and polynomial in at most seven symbolic parameters, one fewer than
  [`GKSLNormalForm`](@ref) because ``\\pi`` takes a parameter slot, and the drive frequency may
  not appear in a coefficient.
- The jump times are named ``u_1, \\dots, u_N``, so no parameter may carry these names.
- The cost grows with the number of jump sequences, which is exponential in `order` for several
  channels.

# Examples

```jldoctest
julia> h = NLevelSpace(:qubit, 2); σ(i, j) = Transition(h, :σ, i, j);

julia> @variables ω::Real t::Real;

julia> H = (1 // 2) * (σ(1, 1) - σ(2, 2)) + cos(ω * t) * (σ(1, 2) + σ(2, 1));

julia> Φ = period_map(H, ω, t, 2; channels=(jump(σ(1, 2), 1 // 5),))
PeriodMap of order 2

julia> iszero(trace_defect(Φ))
false
```
"""
function period_map(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  order::Int;
  channels::LiouvillianChannelCollection=(),
  algorithm::Union{BlochFeshbach,HoriDeprit}=HoriDeprit(),
)
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  hamiltonian_harmonics = harmonics(qadd(H), wd, t)
  require_hermitian_drive(hamiltonian_harmonics)
  provenance = microscopic_provenance(channels, wd, t)
  generator = harmonics(liouvillian_from_provenance(H, provenance), wd, t)
  lowering, coefficients, twopi, times = period_map_lowering(generator, order)
  identity = algebra_identity(lowering.algebra)
  zero_element = zero(identity)

  lowered = lowered_channels(lowering, provenance, wd, t)
  Q = no_jump_generator(
    lower_harmonics(lowering, hamiltonian_harmonics), lowered, zero_element
  )
  frame = no_jump_frame(algorithm, Q, order)
  drive_scales = [convert(SQA.CNum, inverse_drive_power(wd, n)) for n in 0:(2 * order)]
  lift = CoefficientLift(lowering, coefficients)
  records = JumpRecords(
    [transported_jump(frame, channel.operator, order) for channel in lowered],
    effective_powers(frame, order),
    graded_value(frame.kick, zero_element),
    graded_value(frame.inverse, zero_element),
    [channel.rate for channel in lowered],
    twopi,
    SimplexIntegrals(inv(twopi)),
    lift,
    drive_scales,
    times,
    Dict{TimeMonomial,SQA.CNum}(),
  )

  amplitudes = JumpAmplitude[]
  series = [zero(sandwich(identity, identity)) for _ in 0:(2 * order)]
  for q in 0:order
    q > 0 && isempty(lowered) && break
    R = order - q
    kernel = drift(kernel_identity(identity, q), records.powers, twopi, 1, 0, R)
    jump_records!(series, amplitudes, records, kernel, Int[], q, R)
  end
  defect = graded_trace_defect(series, identity, order)

  return PeriodMap(
    lift_graded_superoperator(lift, series, drive_scales),
    amplitudes,
    lift_scaled_sum(lift, zip(defect, drive_scales)),
    order,
  )
end
