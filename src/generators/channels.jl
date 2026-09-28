const LiouvillianChannelCollection = Union{Tuple,AbstractVector}

abstract type FloquetProvenance end
struct NoProvenance <: FloquetProvenance end

struct NonnegativeRateAssumption
  rate::SQA.CNum
end

struct MicroscopicProvenance <: FloquetProvenance
  collapse_operators::Vector{SQA.QAdd}
  jump_operators::Vector{SQA.QAdd}
  jump_rates::Vector{SQA.CNum}
  rate_assumptions::Vector{NonnegativeRateAssumption}
  frame_seeds::Vector{SQA.QAdd}
end

abstract type LiouvillianChannel end

struct CollapseChannel{O<:SQA.QField} <: LiouvillianChannel
  operator::O
end

struct RateWeightedJump{O<:SQA.QField} <: LiouvillianChannel
  operator::O
  rate::SQA.CNum
  assumption::NonnegativeRateAssumption
end

const MAX_DISPLAYED_CHANNELS = 4

function Base.show(io::IO, channel::CollapseChannel)
  print(io, "collapse(")
  show(io, channel.operator)
  return print(io, ")")
end

function Base.show(io::IO, channel::RateWeightedJump)
  print(io, "jump(")
  show(io, channel.operator)
  print(io, ", ")
  show(io, channel.rate)
  return print(io, ")")
end

function show_channel_term(io::IO, channel::CollapseChannel)
  print(io, "𝒟[")
  show(IOContext(io, :compact => true), channel.operator)
  return print(io, "]")
end

function show_channel_term(io::IO, channel::RateWeightedJump)
  print(io, "(")
  show(IOContext(io, :compact => true), channel.rate)
  print(io, ")𝒟[")
  show(IOContext(io, :compact => true), channel.operator)
  return print(io, "]")
end

function Base.show(io::IO, ::MIME"text/plain", channel::CollapseChannel)
  return show_channel_term(io, channel)
end
function Base.show(io::IO, ::MIME"text/plain", channel::RateWeightedJump)
  return show_channel_term(io, channel)
end

function latex_fragment(x)
  body = strip(sprint(show, MIME"text/latex"(), x))
  for (opening, closing) in
      (("\$\$", "\$\$"), ("\$", "\$"), (raw"\begin{equation}", raw"\end{equation}"))
    if startswith(body, opening) && endswith(body, closing)
      body = strip(chopsuffix(chopprefix(body, opening), closing))
    end
  end
  return body
end

function show_channel_latex(io::IO, channel::CollapseChannel)
  print(io, raw"\mathcal{D}\!\left[")
  print(io, latex_fragment(channel.operator))
  return print(io, raw"\right]")
end

function show_channel_latex(io::IO, channel::RateWeightedJump)
  print(io, raw"\left(")
  print(io, latex_fragment(channel.rate))
  print(io, raw"\right)\mathcal{D}\!\left[")
  print(io, latex_fragment(channel.operator))
  return print(io, raw"\right]")
end

function Base.show(io::IO, ::MIME"text/latex", channel::LiouvillianChannel)
  print(io, raw"\[")
  show_channel_latex(io, channel)
  return print(io, raw"\]")
end

channel_collection_label(::Type{<:CollapseChannel}) = "collapse channel"
channel_collection_label(::Type{<:RateWeightedJump}) = "rate-weighted jump channel"
channel_collection_label(::Type{<:LiouvillianChannel}) = "channel"

function displayed_channel_indices(io::IO, count::Int)
  count == 0 && return Int[]
  get(io, :limit, false) || return collect(1:count)

  rows, _ = displaysize(io)
  maximum_items = min(MAX_DISPLAYED_CHANNELS, max(1, rows - 3))
  count <= maximum_items && return collect(1:count)

  leading = cld(maximum_items, 2)
  trailing = maximum_items - leading
  indices = collect(1:leading)
  trailing > 0 && append!(indices, (count - trailing + 1):count)
  return indices
end

function show_omitted_channels(io::IO, count::Int)
  print(io, "\n  ⋮ ", count, count == 1 ? " channel omitted" : " channels omitted")
  return nothing
end

function Base.show(
  io::IO, ::MIME"text/plain", channels::Vector{T}
) where {T<:LiouvillianChannel}
  count = length(channels)
  label = channel_collection_label(T)
  print(io, count, " ", label, count == 1 ? "" : "s")
  isempty(channels) && return nothing
  print(io, ":")

  previous = 0
  for index in displayed_channel_indices(io, count)
    omitted = index - previous - 1
    omitted > 0 && show_omitted_channels(io, omitted)
    print(io, "\n  ", index, ": ")
    show_channel_term(io, channels[index])
    previous = index
  end
  trailing = count - previous
  trailing > 0 && show_omitted_channels(io, trailing)
  return nothing
end

function show_omitted_channels_latex(io::IO, count::Int)
  print(
    io,
    raw"\\&\quad\vdots\quad\text{(",
    count,
    count == 1 ? " channel omitted" : " channels omitted",
    raw")}",
  )
  return nothing
end

function Base.show(io::IO, ::MIME"text/latex", channels::Vector{<:LiouvillianChannel})
  isempty(channels) && return print(io, raw"\[\mathcal{L}_{\mathrm{diss}} = 0\]")

  print(io, raw"\[\begin{aligned}")
  previous = 0
  first_term = true
  for index in displayed_channel_indices(io, length(channels))
    omitted = index - previous - 1
    omitted > 0 && show_omitted_channels_latex(io, omitted)
    first_term ? print(io, raw"\mathcal{L}_{\mathrm{diss}} &= ") : print(io, raw"\\&+ ")
    show_channel_latex(io, channels[index])
    first_term = false
    previous = index
  end
  trailing = length(channels) - previous
  trailing > 0 && show_omitted_channels_latex(io, trailing)
  return print(io, raw"\end{aligned}\]")
end

function append_frame_seeds!(
  seeds::Vector{SQA.QAdd}, operator::SQA.QAdd, wd::Symbolics.Num, t::Symbolics.Num
)
  lowered = harmonics(operator, wd, t)
  for harmonic in sort!(collect(keys(lowered)); by=label -> (abs(label), label))
    push!(seeds, lowered[harmonic])
  end
  return seeds
end

function microscopic_provenance(
  channels::LiouvillianChannelCollection, wd::Symbolics.Num, t::Symbolics.Num
)
  collapse_operators = SQA.QAdd[]
  jump_operators = SQA.QAdd[]
  jump_rates = SQA.CNum[]
  rate_assumptions = NonnegativeRateAssumption[]
  frame_seeds = SQA.QAdd[]

  for channel in channels
    if channel isa CollapseChannel
      push!(collapse_operators, qadd(channel.operator))
      append_frame_seeds!(frame_seeds, last(collapse_operators), wd, t)
    elseif channel isa RateWeightedJump
      push!(jump_operators, qadd(channel.operator))
      push!(jump_rates, channel.rate)
      push!(rate_assumptions, channel.assumption)
      append_frame_seeds!(frame_seeds, last(jump_operators), wd, t)
    else
      throw(
        ArgumentError("channels must contain only `collapse(...)` and `jump(...)` values")
      )
    end
  end

  return MicroscopicProvenance(
    collapse_operators, jump_operators, jump_rates, rate_assumptions, frame_seeds
  )
end

function liouvillian_from_provenance(H::SQA.QField, provenance::MicroscopicProvenance)
  generator = hamiltonian_action(H)
  for operator in provenance.collapse_operators
    generator = generator + dissipator(operator)
  end
  for i in eachindex(provenance.jump_operators)
    generator =
      generator + provenance.jump_rates[i] * dissipator(provenance.jump_operators[i])
  end
  return generator
end

@inline channel_liouvillian(channel::CollapseChannel) = dissipator(channel.operator)
@inline channel_liouvillian(channel::RateWeightedJump) =
  channel.rate * dissipator(channel.operator)
"""
    liouvillian(H::QField; channels=()) -> Liouvillian

Construct
``ρ ↦ -i[H, ρ] + Σₐ D[Cₐ](ρ) + Σᵦ γᵦD[Jᵦ](ρ)``.

`channels` is a tuple or vector of [`collapse`](@ref) and [`jump`](@ref) values. See
[`jump`](@ref) for the physical-rate requirements on `γᵦ`.

# Arguments

- `H`: Symbolic Hamiltonian expression.
- `channels`: Tuple or vector of channel values to add to the coherent action.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables γ::Real;

julia> L = liouvillian(a' * a; channels=(jump(a, γ),));

julia> L == hamiltonian_action(a' * a) + γ * dissipator(a)
true
```
"""
function liouvillian(H::SQA.QField; channels::LiouvillianChannelCollection=())
  generator = hamiltonian_action(H)
  for channel in channels
    channel isa LiouvillianChannel || throw(
      ArgumentError("channels must contain only `collapse(...)` and `jump(...)` values")
    )
    generator = generator + channel_liouvillian(channel)
  end
  return generator
end

"""
    collapse(operator::QField) -> CollapseChannel

Create a channel from a complete collapse operator `L`. When passed to a
[`liouvillian`](@ref), it contributes the [`dissipator`](@ref) of `operator` with unit weight:
``D[L](ρ) = LρL† - (L†Lρ + ρL†L)/2``.
Here, `L` denotes `operator`; any amplitude or phase belonging to the collapse process is included
in it.

See also [`jump`](@ref), [`liouvillian`](@ref).
"""
function collapse(operator::SQA.QField)
  return CollapseChannel(operator)
end

function known_numeric_jump_rate(rate::LiouvillianScalar)
  if rate isa Symbolics.Num
    value = Symbolics.value(rate)
    return value isa Real ? value : nothing
  elseif rate isa Real
    return rate
  elseif rate isa Complex && iszero(imag(rate))
    real_rate = real(rate)
    if real_rate isa Symbolics.Num
      value = Symbolics.value(real_rate)
      return value isa Real ? value : nothing
    elseif real_rate isa Real
      return real_rate
    end
  end
  return nothing
end

function validated_jump_rate(rate::LiouvillianScalar)
  coefficient = convert(SQA.CNum, rate)
  real_valued =
    rate isa Real ||
    (rate isa Complex && iszero(imag(rate))) ||
    coefficient == conj(coefficient)
  real_valued || throw(ArgumentError("jump rate must be provably real; got `$rate`"))

  numeric_rate = known_numeric_jump_rate(rate)
  numeric_rate !== nothing &&
    numeric_rate < 0 &&
    throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))

  return coefficient
end

"""
    jump(operator::QField, rate) -> RateWeightedJump

Create a physical rate-weighted channel ``γ D[J]`` from a bare jump operator `J` and rate `γ`.
The rate may be time dependent and periodic, but it must be provably real. Negative numeric rates
are rejected. A real symbolic expression is accepted under the assumption that the expression as
a whole is nonnegative; no sign analysis of its factors is performed.

Thus `jump(J, γ₀ + γ₁*cos(ω*t))` is supported for real symbolic parameters. A genuinely complex
coefficient such as `expim(ω*t)` is not a physical rate. For signed or complex algebraic
coefficients, use `c * dissipator(J)` instead.

Use [`collapse`](@ref) when the complete channel amplitude is already folded into the operator.

See also [`collapse`](@ref), [`dissipator`](@ref), [`liouvillian`](@ref).
"""
function jump(operator::SQA.QField, rate::LiouvillianScalar)
  coefficient = validated_jump_rate(rate)
  return RateWeightedJump(operator, coefficient, NonnegativeRateAssumption(coefficient))
end
