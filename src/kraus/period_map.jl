"""
    JumpAmplitude

One continuous Kraus operator of a [`PeriodMapExpansion`](@ref): the retained amplitude ``K(u)`` of the
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

The exact one-period map of a periodically driven Lindblad system, from ``t = 0`` to
``t = T = 2\\pi/\\omega_d``, returned by [`period_map`](@ref). It holds the Hamiltonian and the
jump channels and is not expanded; [`floquet_expansion`](@ref) expands it in powers of
``1/\\omega_d`` into a [`PeriodMapExpansion`](@ref).
"""
struct PeriodMap
  hamiltonian::PeriodicGenerator{SQA.QAdd}
  operators::Vector{PeriodicGenerator{SQA.QAdd}}
  rates::Vector{SQA.CNum}
  generator::PeriodicGenerator{Liouvillian}
end

function Base.show(io::IO, ::MIME"text/plain", Φ::PeriodMap)
  n = length(Φ.operators)
  return print(io, "PeriodMap with ", n, n == 1 ? " channel" : " channels")
end
Base.show(io::IO, Φ::PeriodMap) = show(io, MIME"text/plain"(), Φ)

"""
    PeriodMapExpansion

Result of [`floquet_expansion`](@ref) applied to a [`PeriodMap`](@ref): a completely positive map
``\\Phi`` that agrees with the exact evolution from ``t = 0`` to ``t = T = 2\\pi/\\omega`` through
``\\omega^{-N}``, where ``N`` is [`order`](@ref). Read it with [`superoperator`](@ref),
[`jump_amplitudes`](@ref) and [`trace_defect`](@ref) rather than by field access.
"""
struct PeriodMapExpansion{G<:Gauge}
  superoperator::Liouvillian
  amplitudes::Vector{JumpAmplitude}
  trace_defect::SQA.QAdd
  gauge::G
  order::Int
end

function Base.show(io::IO, ::MIME"text/plain", Φ::PeriodMapExpansion{G}) where {G}
  return print(io, "PeriodMapExpansion{", G, "} of order ", Φ.order)
end
Base.show(io::IO, Φ::PeriodMapExpansion) = show(io, MIME"text/plain"(), Φ)

"""
    order(Φ::PeriodMapExpansion) -> Int

Return the inverse-frequency order ``N`` through which `Φ` agrees with the exact period map.
"""
order(Φ::PeriodMapExpansion) = Φ.order

"""
    superoperator(Φ::PeriodMapExpansion) -> Liouvillian

Return the period map as the linear map ``\\rho \\mapsto \\sum_j c_j A_j \\rho B_j``. Its
coefficients are exact polynomials in ``\\pi``, the parameters, and ``1/\\omega``. The map is
completely positive for every drive frequency and nonnegative jump rates, not only through the
retained order, because it is a sum over jump records of squared jump amplitudes; see
[Floquet–Kraus period map](@ref floquet-kraus-theory).

`SQA.simplify` evaluates ``\\pi`` in floating point, so simplify a coefficient only after
substituting values, or substitute with `fold=Val(true)`.
"""
superoperator(Φ::PeriodMapExpansion) = Φ.superoperator

"""
    trace_defect(Φ::PeriodMapExpansion) -> QAdd

Return ``\\Phi^\\dagger(1) - 1``, which measures the failure of trace preservation. It starts at
``\\omega^{-N-1}``, where ``N`` is [`order`](@ref); the construction checks exactly that every
lower order vanishes. The map ``\\rho \\mapsto \\Phi(S\\rho S)`` with ``S = (1 + D)^{-1/2}``, ``D``
the trace defect, is exactly trace preserving and still agrees with the exact period map through
``\\omega^{-N}``, wherever ``1 + D`` is positive definite.
"""
trace_defect(Φ::PeriodMapExpansion) = Φ.trace_defect

"""
    jump_amplitudes(Φ::PeriodMapExpansion) -> Vector{JumpAmplitude}

Return the period map in continuous Kraus form,
``\\Phi(\\rho) = \\sum_K w_K \\int K(u) \\rho K(u)^\\dagger \\, du``, with one
[`JumpAmplitude`](@ref) per sequence of jump channels. The integral over the jump times reproduces
[`superoperator`](@ref) exactly. A finite Kraus form follows at substituted values from a
factorization of the Gram matrix of the monomials in ``u``.
"""
jump_amplitudes(Φ::PeriodMapExpansion) = Φ.amplitudes

struct LoweredChannel{X,F}
  operator::Harmonics{X}
  rate::F
end

function lower_harmonics(lowering::SQALowering, P::PeriodicGenerator{SQA.QAdd})
  return Dict(k => lower_qadd(lowering, P[k]) for k in keys(P))
end

function check_time_independent_rate(rate::SQA.CNum, wd::Symbolics.Num, t::Symbolics.Num)
  found = Symbolics.get_variables(SQA.to_num(rate))
  any(v -> isequal(v, Symbolics.unwrap(t)) || isequal(v, Symbolics.unwrap(wd)), found) &&
    throw(
      ArgumentError(
        "the period map needs time-independent jump rates, found $rate; put the time " *
        "dependence into the jump operator",
      ),
    )
  return rate
end

function lowered_channels(lowering::SQALowering{T}, Φ::PeriodMap) where {T}
  return LoweredChannel{AlgebraOperator{T,real(T)},T}[
    LoweredChannel(lower_harmonics(lowering, operator), lower_coefficient(lowering, rate))
    for (operator, rate) in zip(Φ.operators, Φ.rates)
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
  filter!(entry -> !iszero(last(entry)), Q)
  haskey(Q, 0) || (Q[0] = zero_element)
  return Q
end

function period_map_scalars(
  generator::PeriodicGenerator{Liouvillian}, ::Tuple{}, ::Type{T}
) where {T}
  return SQALowering{T}(generator), ExactCoefficients()
end

function period_map_scalars(
  generator::PeriodicGenerator{Liouvillian},
  parameters::Tuple{Symbolics.Num,Vararg{Symbolics.Num}},
  ::Type,
)
  T = Complex{ParameterScalar}
  symbols = Symbolics.Num[p for p in parameters]
  values = T[parameter_variable(i) for i in eachindex(symbols)]
  lowering = SQALowering{T}(generator, symbols, values)
  return lowering, ParameterCoefficients(symbols, ones(Int, length(symbols)))
end

function period_map_lowering(
  generator::PeriodicGenerator{Liouvillian}, order::Int, ::Type{T}=Complex{Rational{BigInt}}
) where {T}
  parameters = checked_native_parameters(generator, MAXIMUM_PARAMETERS, "the period map")
  times = [Symbolics.variable(:u, j; T=Real) for j in 1:order]
  any(u -> any(p -> isequal(p, u), parameters), times) && throw(
    ArgumentError(
      "the jump times are named u₁, …, u_$order; rename the parameters that use these names",
    ),
  )
  return period_map_scalars(generator, Tuple(parameters), T)..., times
end

struct JumpRecords{X,F,L<:CoefficientLift}
  jumps::Vector{Vector{Harmonics{X}}}
  powers::Vector{Vector{X}}
  left::Vector{X}
  right::Vector{X}
  rates::Vector{F}
  integrals::SimplexIntegrals{F}
  lift::L
  drive_scales::Vector{SQA.CNum}
  period_scales::Vector{SQA.CNum}
  times::Vector{Symbolics.Num}
  time_values::Dict{TimeMonomial,SQA.CNum}
end

function record_rate(word::Vector{Int}, records::JumpRecords{X,F}) where {X,F}
  rate = one(F)
  for channel in word
    rate *= records.rates[channel]
  end
  return rate
end

function period_scale(records::JumpRecords, k::Int)
  scales = records.period_scales
  while length(scales) <= k
    push!(scales, last(scales) * convert(SQA.CNum, 2 * Symbolics.Num(π)))
  end
  return scales[k + 1]
end

const GradedSuperoperators{S} = Dict{Tuple{Int,Int},S}

struct RecordTargets{S}
  zero::S
  diagonal::GradedSuperoperators{S}
  upper::GradedSuperoperators{S}
  amplitudes::Vector{JumpAmplitude}
end

function RecordTargets(identity::AlgebraOperator)
  zero_superoperator = zero(sandwich(identity, identity))
  return RecordTargets(
    zero_superoperator,
    GradedSuperoperators{typeof(zero_superoperator)}(),
    GradedSuperoperators{typeof(zero_superoperator)}(),
    JumpAmplitude[],
  )
end

function graded_target!(targets::RecordTargets, graded, key)
  return get!(() -> zero(targets.zero), graded, key)
end

function upper_gram_row(
  entries::Vector, adjoints::Vector{X}, i::Int, rate::T, integrals
) where {T,X}
  (_, ℓ, μ) = first(entries[i])
  weighted = Dict{Tuple{Int,Int},Dict{Monomial,T}}()
  for j in (i + 1):length(entries)
    s, m, ν = first(entries[j])
    for (k, g) in gram_entry(integrals, μ, ν)
      terms = get!(Dict{Monomial,T}, weighted, (s, ℓ + m + k))
      add_scaled!(terms, rate * g, adjoints[j])
    end
  end
  return weighted
end

function record_map!(
  targets::RecordTargets,
  q::Int,
  kernel::TimeKernel{AlgebraOperator{T,R}},
  rate::T,
  integrals,
) where {T,R}
  entries = collect(kernel.terms)
  adjoints = [adjoint(B) for (_, B) in entries]
  for (i, ((r, ℓ, μ), A)) in pairs(entries)
    for (k, g) in gram_entry(integrals, μ, μ)
      target = graded_target!(targets, targets.diagonal, (q + 2r, 2ℓ + k))
      add_sandwich!(target, A, adjoints[i], rate * g)
    end
    for ((s, k), terms) in upper_gram_row(entries, adjoints, i, rate, integrals)
      target = graded_target!(targets, targets.upper, (q + r + s, k))
      add_sandwich!(target, A, AlgebraOperator{T,R}(A.algebra, terms))
    end
  end
  return targets
end

function time_monomial((phases, powers)::TimeMonomial, times::Vector{Symbolics.Num})
  value = Complex{Symbolics.Num}(1)
  for (u, k, p) in zip(times, phases, powers)
    angle = 2 * Symbolics.Num(π) * k * u
    value *= u^p * Complex(cos(angle), -sin(angle))
  end
  return convert(SQA.CNum, value)
end

function time_value(records::JumpRecords, monomial::TimeMonomial)
  return get!(() -> time_monomial(monomial, records.times), records.time_values, monomial)
end

function lift_amplitude(records::JumpRecords, kernel::TimeKernel, rate, q::Int)
  weight =
    lift_coefficient(records.lift.coefficient, rate) *
    records.drive_scales[q + 1] *
    period_scale(records, q)
  pieces = (
    (
      X,
      records.drive_scales[r + 1] *
      period_scale(records, ℓ + sum(last(monomial); init=0)) *
      time_value(records, monomial),
    ) for ((r, ℓ, monomial), X) in kernel.terms
  )
  return JumpAmplitude(weight, lift_scaled_sum(records.lift, pieces), records.times[1:q])
end

function record!(targets::RecordTargets, records::JumpRecords, kernel, word, R::Int)
  amplitude = conjugate_kernel(kernel, records.left, records.right, R)
  isempty(amplitude.terms) && return targets
  rate = record_rate(word, records)
  record_map!(targets, length(word), amplitude, rate, records.integrals)
  push!(targets.amplitudes, lift_amplitude(records, amplitude, rate, length(word)))
  return targets
end

function jump_records!(
  targets::RecordTargets, records::JumpRecords, kernel, word::Vector{Int}, q::Int, R::Int
)
  isempty(kernel.terms) && return targets
  j = length(word)
  j == q && return record!(targets, records, kernel, word, R)
  for (channel, jump) in pairs(records.jumps)
    next = insert_jump(kernel, jump, j + 1, R)
    push!(word, channel)
    jump_records!(
      targets, records, drift(next, records.powers, j + 2, j + 1, R), word, q, R
    )
    pop!(word)
  end
  return targets
end

function dual_identity(S::AlgebraSuperoperator{T,R}) where {T,R}
  terms = Dict{Monomial,T}()
  for ((left, right), c) in S.terms,
    (monomial, d) in monomial_product(S.algebra, right, left)

    accumulate!(terms, monomial, c * d)
  end
  return AlgebraOperator{T,R}(S.algebra, terms)
end

function graded_trace_defect(series::GradedSuperoperators, identity, N::Int)
  defect = Dict(key => dual_identity(S) for (key, S) in series)
  defect[(0, 0)] -= identity
  all(iszero(D) for ((n, _), D) in defect if n <= N) ||
    error("the period map failed its trace-preservation certificate")
  return defect
end

function lift_graded(records::JumpRecords, graded::Dict)
  return (
    (X, records.drive_scales[n + 1] * period_scale(records, k)) for
    ((n, k), X) in graded if !iszero(X)
  )
end

function lift_graded_superoperator(
  lift::CoefficientLift, records::JumpRecords, series::GradedSuperoperators
)
  result = Liouvillian(LiouvillianTerms())
  for (S, scale) in lift_graded(records, series)
    lift_superoperator!(result, lift, S, scale)
  end
  return result
end

"""
    period_map(H::Union{QSym,QAdd}, ωd, t; channels=()) -> PeriodMap

Return the one-period map of the periodically driven Lindblad system

```math
\\dot\\rho = -i[H(t), \\rho] + \\sum_a \\gamma_a \\mathcal{D}[L_a(t)]\\rho
```

from ``t = 0`` to ``t = T = 2\\pi/\\omega_d``. The map is exact and not expanded; expand it in
powers of ``1/\\omega_d`` with [`floquet_expansion`](@ref). It keeps the jump operators, which the
Liouvillian alone does not fix, because the expansion is built from quantum-jump records.

# Arguments

- `H`: Symbolic time-dependent Hamiltonian, Hermitian at every time.
- `ωd`: Symbolic drive frequency.
- `t`: Symbolic time variable.
- `channels`: Tuple or vector of [`collapse`](@ref) and [`jump`](@ref) values.

Jump rates must be time independent; a periodic rate belongs in the operator.

# Examples

```jldoctest
julia> h = NLevelSpace(:qubit, 2); σ(i, j) = Transition(h, :σ, i, j);

julia> @variables ω::Real t::Real;

julia> H = (1 // 2) * (σ(1, 1) - σ(2, 2)) + cos(ω * t) * (σ(1, 2) + σ(2, 1));

julia> Φ = period_map(H, ω, t; channels=(jump(σ(1, 2), 1 // 5),))
PeriodMap with 1 channel
```
"""
function period_map(
  H::OperatorExpression,
  wd::Symbolics.Num,
  t::Symbolics.Num;
  channels::LiouvillianChannelCollection=(),
)
  hamiltonian = harmonics(qadd(H), wd, t)
  require_hermitian_drive(hamiltonian)
  provenance = microscopic_provenance(channels, wd, t)
  operators = [provenance.collapse_operators; provenance.jump_operators]
  rates = SQA.CNum[
    fill(convert(SQA.CNum, 1), length(provenance.collapse_operators))
    [check_time_independent_rate(rate, wd, t) for rate in provenance.jump_rates]
  ]
  return PeriodMap(
    hamiltonian,
    [harmonics(operator, wd, t) for operator in operators],
    rates,
    harmonics(liouvillian_from_provenance(H, provenance), wd, t),
  )
end

"""
    floquet_expansion(Φ::PeriodMap, gauge::VanVleck, order) -> PeriodMapExpansion

Expand the one-period map `Φ` in powers of ``1/\\omega_d`` as a completely positive map. The
result agrees with the exact period map through ``\\omega_d^{-\\text{order}}``, so its logarithm
divided by ``T`` agrees with the stroboscopic effective generator through
``\\omega_d^{1-\\text{order}}``, as for the generator expansion with the same `order`.

The map is never assumed to have a time-independent GKSL logarithm. It is the sum over
quantum-jump records of the period of squared jump amplitudes, integrated exactly over the jump
times. Every amplitude is retained through the requested order and squared without truncation, so
the map is completely positive for every drive frequency, and its [`trace_defect`](@ref) starts
beyond the retained order. A symbolic rate enters as the nonnegative quantity that [`jump`](@ref)
assumes it to be, and complete positivity holds only where that assumption does.

`gauge` fixes the Floquet frame of the no-jump evolution, and its algorithm,
[`HoriDeprit`](@ref) or [`BlochFeshbach`](@ref), computes that frame; both algorithms give the same
map. Only [`VanVleck`](@ref) is supported. The construction is described in
[Floquet–Kraus period map](@ref floquet-kraus-theory).

# Restrictions

- Coefficients must be exact and polynomial in at most eight symbolic parameters, as for
  [`GKSLNormalForm`](@ref), and the drive frequency may not appear in a coefficient.
- The jump times are named ``u_1, \\dots, u_N``, so no parameter may carry these names.
- The cost grows with the number of jump sequences, which is exponential in `order` for several
  channels.

# Examples

```jldoctest
julia> h = NLevelSpace(:qubit, 2); σ(i, j) = Transition(h, :σ, i, j);

julia> @variables ω::Real t::Real;

julia> H = (1 // 2) * (σ(1, 1) - σ(2, 2)) + cos(ω * t) * (σ(1, 2) + σ(2, 1));

julia> Φ = period_map(H, ω, t; channels=(jump(σ(1, 2), 1 // 5),));

julia> x = floquet_expansion(Φ, VanVleck(), 2)
PeriodMapExpansion{VanVleck{HoriDeprit}} of order 2

julia> iszero(trace_defect(x))
false
```
"""
function floquet_expansion(Φ::PeriodMap, gauge::VanVleck, order::Int)
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  try
    return expand_period_map(
      Φ, gauge, order, period_map_lowering(Φ.generator, order, Complex{Rational{Int128}})...
    )
  catch error
    error isa Union{OverflowError,InexactError} || rethrow()
    return expand_period_map(Φ, gauge, order, period_map_lowering(Φ.generator, order)...)
  end
end

function expand_period_map(
  Φ::PeriodMap,
  gauge::VanVleck,
  order::Int,
  lowering::SQALowering{T},
  coefficients,
  times::Vector{Symbolics.Num},
) where {T}
  wd = Φ.generator.wd
  identity = algebra_identity(lowering.algebra)
  zero_element = zero(identity)

  lowered = lowered_channels(lowering, Φ)
  Q = no_jump_generator(lower_harmonics(lowering, Φ.hamiltonian), lowered, zero_element)
  frame = no_jump_frame(gauge.algorithm, Q, order)
  drive_scales = [convert(SQA.CNum, inverse_drive_power(wd, n)) for n in 0:(2 * order)]
  lift = CoefficientLift(lowering, coefficients)
  records = JumpRecords(
    [transported_jump(frame, channel.operator, order) for channel in lowered],
    effective_powers(frame, order),
    endpoint_value(frame.kick, zero_element),
    endpoint_value(frame.inverse, zero_element),
    T[channel.rate for channel in lowered],
    SimplexIntegrals{T}(),
    lift,
    drive_scales,
    [convert(SQA.CNum, 1)],
    times,
    Dict{TimeMonomial,SQA.CNum}(),
  )

  targets = RecordTargets(identity)
  for q in 0:order
    q > 0 && isempty(lowered) && break
    R = order - q
    kernel = drift(kernel_identity(identity, q), records.powers, 1, 0, R)
    jump_records!(targets, records, kernel, Int[], q, R)
  end
  series = GradedSuperoperators{typeof(targets.zero)}()
  for key in union(keys(targets.diagonal), keys(targets.upper))
    D = get(targets.diagonal, key, targets.zero)
    U = get(targets.upper, key, targets.zero)
    series[key] = add_adjoint_mirror!(D + U, U)
  end
  defect = graded_trace_defect(series, identity, order)

  return PeriodMapExpansion(
    lift_graded_superoperator(lift, records, series),
    targets.amplitudes,
    lift_scaled_sum(lift, lift_graded(records, defect)),
    gauge,
    order,
  )
end

function floquet_expansion(::PeriodMap, gauge::Gauge, ::Int)
  return throw(
    ArgumentError(
      "a period map is expanded in the van Vleck gauge only, found $(typeof(gauge)); " *
      "pass VanVleck()",
    ),
  )
end
