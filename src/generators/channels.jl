abstract type FloquetProvenance end
struct NoProvenance <: FloquetProvenance end

struct NonnegativeRateAssumption
  rate::SQA.CNum
end

struct MicroscopicProvenance <: FloquetProvenance
  jump_operators::Vector{SQA.QAdd}
  jump_rates::Vector{SQA.CNum}
  rate_assumptions::Vector{NonnegativeRateAssumption}
  frame_seeds::Vector{SQA.QAdd}
end

function append_frame_seeds!(seeds::Vector{SQA.QAdd}, lowered::PeriodicGenerator{SQA.QAdd})
  for harmonic in sort!(collect(keys(lowered)); by=label -> (abs(label), label))
    push!(seeds, lowered[harmonic])
  end
  return seeds
end

function known_negative_jump_rate(coefficient::SQA.CNum)
  SQA.has_coefficient_terms(coefficient) || return false
  terms = SQA.coefficient_terms(coefficient)
  return !isempty(terms) && all(terms) do term
    return isempty(term.powers) && isreal(term.scalar) && real(term.scalar) < 0
  end
end

function validated_jump_rate(coefficient::SQA.CNum)
  coefficient == conj(coefficient) ||
    throw(ArgumentError("channel rate must be provably real; got `$coefficient`"))
  known_negative_jump_rate(coefficient) &&
    throw(ArgumentError("channel rate must be nonnegative; got `$coefficient`"))
  return coefficient
end

function microscopic_jump(channel::SQA.Channel{SQA.QAdd})
  isempty(SQA.scope(channel)) || throw(
    ArgumentError(
      "a channel summed over a symbolic index has no microscopic rate per site; " *
      "expand `SuperOperator(L)` instead",
    ),
  )
  operator = SQA.jump(channel)
  isequal(operator, SQA.partner(channel)) || throw(
    ArgumentError(
      "a correlated channel 𝒟[J, K] has no single nonnegative rate; " *
      "expand `SuperOperator(L)` instead",
    ),
  )
  return operator
end

function microscopic_provenance(
  L::SQA.Lindbladian{SQA.QAdd}, wd::Symbolics.Num, t::Symbolics.Num
)
  jump_operators = SQA.QAdd[]
  jump_rates = SQA.CNum[]
  rate_assumptions = NonnegativeRateAssumption[]
  frame_seeds = SQA.QAdd[]

  for (channel, rate) in SQA.channels(L)
    operator = microscopic_jump(channel)
    coefficient = validated_jump_rate(rate)
    push!(jump_operators, operator)
    push!(jump_rates, coefficient)
    push!(rate_assumptions, NonnegativeRateAssumption(coefficient))
    append_frame_seeds!(frame_seeds, harmonics(operator, wd, t))
  end

  return MicroscopicProvenance(jump_operators, jump_rates, rate_assumptions, frame_seeds)
end

function gksl_lindbladian(H::SQA.QAdd, jumps)
  generator = SQA.hamiltonian_action(H)
  for (operator, rate) in jumps
    generator = generator + rate * SQA.dissipator(operator)
  end
  return generator
end
