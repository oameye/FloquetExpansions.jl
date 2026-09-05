"""
    FloquetExpansion

Result of [`floquet_expansion`](@ref). Read it with [`effective_generator`](@ref),
[`effective_component`](@ref), and [`micromotion`](@ref) rather than by field access.

Stored expansion coefficients do not include inverse powers of the drive frequency; the public
accessors restore that scaling. Positive completion may change [`effective_generator`](@ref), but
never the retained [`effective_component`](@ref) or [`micromotion`](@ref).

See also [`floquet_expansion`](@ref), [`positive_completion`](@ref).
"""
struct FloquetExpansion{G<:Gauge,P<:PeriodicGenerator,E,C<:Completion,R<:FloquetProvenance}
  generator::P
  kick_components::Vector{P}
  effective_components::Vector{E}
  gauge::G
  order::Int
  completion::C
  provenance::R
end

function Base.getproperty(expansion::FloquetExpansion, name::Symbol)
  if name === :provenance
    throw(ArgumentError("FloquetExpansion field :provenance is private"))
  end
  if name === :kick_derivative_components ||
    name === :dressed_generator ||
    name === :dressed_kick_derivative
    throw(
      ArgumentError(
        "FloquetExpansion field :$(name) is private; use `effective_generator`, `effective_component`, and `micromotion`",
      ),
    )
  end
  return getfield(expansion, name)
end

function Base.propertynames(
  ::FloquetExpansion{G,P,E,C,R}, private::Bool=false
) where {G,P,E,C,R}
  names = (:generator, :kick_components, :effective_components, :gauge, :order, :completion)
  return if private
    (
      names...,
      :provenance,
      :kick_derivative_components,
      :dressed_generator,
      :dressed_kick_derivative,
    )
  else
    names
  end
end

const GeneratorComponent = Union{SQA.QAdd,Liouvillian}

triindex(n::Int, j::Int) = (n * (n + 1)) ÷ 2 + j + 1
lie_transform_phase(::Type{SQA.QAdd}) = im
lie_transform_phase(::Type{Liouvillian}) = -1
function weight_generator(::Type{T}, j::Int) where {T<:GeneratorComponent}
  return lie_transform_phase(T)^j * (1 // factorial(j))
end
function weight_kick_derivative(::Type{T}, j::Int) where {T<:GeneratorComponent}
  return lie_transform_phase(T)^j * (1 // factorial(j + 1))
end

function Base.show(io::IO, ::MIME"text/plain", expansion::FloquetExpansion{G}) where {G}
  return print(io, "FloquetExpansion{", nameof(G), "} of order ", expansion.order)
end

Base.show(io::IO, expansion::FloquetExpansion) = show(io, MIME"text/plain"(), expansion)

function dressed_node(K::Vector{P}, previous, n::Int, j::Int, generator::P) where {P}
  result = zero(generator)
  for k in 1:(n - j + 1)
    result = result + SQA.commutator(K[k], previous(k, j))
  end
  return result
end

function assemble_resolvent(
  dressed_generator::Vector{P}, dressed_kick_derivative::Vector{P}, n::Int, generator::P
) where {P<:PeriodicGenerator}
  T = eltype(typeof(generator))
  result = zero(generator)
  for j in 0:n
    result = result + weight_generator(T, j) * dressed_generator[triindex(n, j)]
  end
  for j in 1:n
    result = result - weight_kick_derivative(T, j) * dressed_kick_derivative[triindex(n, j)]
  end
  return result
end

function require_hermitian_drive(generator::PeriodicGenerator{SQA.QAdd})
  LinearAlgebra.ishermitian(generator) || throw(
    ArgumentError(
      "the drive is not Hermitian: it must satisfy H_{-m} = H_m' (eq:fourierH)"
    ),
  )
  return generator
end

function floquet_expansion_impl(
  generator::P, gauge::G, order::Int, provenance::R
) where {P<:PeriodicGenerator,G<:Gauge,R<:FloquetProvenance}
  order >= 1 || throw(ArgumentError("order must be >= 1"))

  generator isa PeriodicGenerator{SQA.QAdd} && require_hermitian_drive(generator)

  nodes = (order * (order + 1)) ÷ 2
  dressed_generator = [zero(generator) for _ in 1:nodes]
  dressed_kick_derivative = [zero(generator) for _ in 1:nodes]
  generator_type = typeof(generator)
  K = generator_type[]
  Kdot = generator_type[]
  E = typeof(time_average(generator))
  effective = E[]

  for n in 0:(order - 1)
    dressed_generator[triindex(n, 0)] = n == 0 ? generator : zero(generator)

    for j in 1:n
      dressed_generator[triindex(n, j)] = dressed_node(
        K, (k, _) -> dressed_generator[triindex(n - k, j - 1)], n, j, generator
      )
    end

    for j in 1:n
      dressed_kick_derivative[triindex(n, j)] = dressed_node(
        K,
        (k, j_) ->
          j_ == 1 ? Kdot[n - k + 1] : dressed_kick_derivative[triindex(n - k, j_ - 1)],
        n,
        j,
        generator,
      )
    end

    resolvent = SQA.simplify(
      assemble_resolvent(dressed_generator, dressed_kick_derivative, n, generator)
    )
    effective_n = SQA.simplify(time_average(resolvent))
    push!(effective, effective_n)

    if n < order - 1
      next_kick = SQA.simplify(antiderivative(remove_average(resolvent), gauge))
      push!(K, next_kick)
      push!(Kdot, derivative(next_kick))
    end
  end

  return FloquetExpansion(generator, K, effective, gauge, order, Uncompleted(), provenance)
end

"""
    floquet_expansion(generator::PeriodicGenerator, gauge, order)
    floquet_expansion(L::Liouvillian, ωd, t, gauge, order)
    floquet_expansion(H::QField, ωd, t, gauge, order; channels=())

Expand a periodically driven generator into a time-independent effective generator and periodic
micromotion, returning a [`FloquetExpansion`](@ref) with `order` retained inverse-frequency
contributions.

# Arguments

- `generator`: Prepared periodic Hamiltonian or Liouvillian generator.
- `L`: Symbolic time-dependent Liouvillian to decompose using `ωd` and `t`.
- `H`: Symbolic time-dependent Hamiltonian to decompose using `ωd` and `t`.
- `ωd`: Symbolic drive frequency.
- `t`: Symbolic time variable.
- `gauge`: Gauge fixing the micromotion integration constant.
- `order`: Number of retained orders; must be at least one.
- `channels`: Tuple or vector of [`collapse`](@ref) and [`jump`](@ref) values added to `H`.

The `L` and `H` forms decompose the time dependence before applying the expansion. Physical
`channels` supplied through the `H` form are retained as microscopic input for positive completion;
explicitly constructed Liouvillian inputs remain valid and can be completed from their Floquet
data or in a supplied [`DissipativeFrame`](@ref).

Hamiltonian generators are checked for Hermiticity at ingest, including the Hamiltonian part of
the `H; channels=...` form. Liouvillian generators use the common algebra without a Hamiltonian
Hermiticity requirement. An explicit channel vector must be nonempty; omit `channels` or pass `()`
for a Hamiltonian-only expansion.

# Notes

`order = 1` retains the period average; the neglected terms start at
``\\mathcal{O}(\\omega_d^{-\\text{order}})``. The high-frequency series is asymptotic rather than
convergent.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables ω::Real t::Real g::Real;

julia> H = harmonics(ω * (a' * a) + g * cos(ω * t) * (a + a'), ω, t);

julia> vv = floquet_expansion(H, VanVleck(), 1)
FloquetExpansion{VanVleck} of order 1

julia> effective_generator(vv)
ω * a' * a

julia> iszero(micromotion(vv))
true
```

See also [`effective_generator`](@ref), [`effective_component`](@ref), [`micromotion`](@ref),
[`positive_completion`](@ref), and [`harmonics`](@ref).
"""
function floquet_expansion(
  generator::P, gauge::G, order::Int
) where {P<:PeriodicGenerator,G<:Gauge}
  return floquet_expansion_impl(generator, gauge, order, NoProvenance())
end

function floquet_expansion(
  L::Liouvillian, wd::Symbolics.Num, t::Symbolics.Num, gauge::Gauge, order::Int
)
  return floquet_expansion_impl(harmonics(L, wd, t), gauge, order, NoProvenance())
end

function floquet_expansion_channels(
  H::SQA.QField, wd::Symbolics.Num, t::Symbolics.Num, gauge::Gauge, order::Int, ::Tuple{}
)
  return floquet_expansion_impl(harmonics(qadd(H), wd, t), gauge, order, NoProvenance())
end

function floquet_expansion_channels(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  gauge::Gauge,
  order::Int,
  channels::LiouvillianChannelCollection,
)
  channels isa AbstractVector &&
    isempty(channels) &&
    throw(
      ArgumentError(
        "an explicit channel vector must be nonempty; omit `channels` for no channels"
      ),
    )

  require_hermitian_drive(harmonics(qadd(H), wd, t))
  provenance = microscopic_provenance(channels)
  L = liouvillian_from_provenance(H, provenance)
  return floquet_expansion_impl(harmonics(L, wd, t), gauge, order, provenance)
end

function floquet_expansion(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  gauge::Gauge,
  order::Int;
  channels::LiouvillianChannelCollection=(),
)
  return floquet_expansion_channels(H, wd, t, gauge, order, channels)
end

function reattach(component::GeneratorComponent, wd::Symbolics.Num, n::Int)
  return iszero(n) ? component : wd^(-n) * component
end
function reattach(generator::PeriodicGenerator{T}, n::Int) where {T<:GeneratorComponent}
  return iszero(n) ? generator : generator.wd^(-n) * generator
end

"""
    effective_generator(expansion::FloquetExpansion) -> T

Return the effective generator represented by `expansion`. Before positive completion this is

``\\mathcal{G}_\\mathrm{eff}^{[N]} = \\sum_{n<N}
\\omega_d^{-n}\\mathcal{G}_\\mathrm{eff}^{(n)}``.

After positive completion it is the completed generator. Use [`effective_component`](@ref) to
inspect the retained perturbative contributions.

See also [`effective_component`](@ref), [`micromotion`](@ref), [`positive_completion`](@ref).
"""
function effective_generator(
  expansion::FloquetExpansion{G,P,E,Uncompleted,R}
) where {G,P,E,R}
  result = zero(expansion.effective_components[1])
  for n in 0:(expansion.order - 1)
    result =
      result + reattach(expansion.effective_components[n + 1], expansion.generator.wd, n)
  end
  return SQA.simplify(result)::E
end

function effective_generator(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P,E,C<:PositiveCompletion,R}
  return getfield(expansion, :completion).generator
end

"""
    effective_component(expansion::FloquetExpansion, n::Int) -> T

Return the retained order-`n` effective-generator contribution, including the scaling
``\\omega_d^{-n}``. The index must satisfy `0 ≤ n < expansion.order`.

Positive completion leaves this contribution unchanged.

See also [`effective_generator`](@ref), [`micromotion`](@ref).
"""
function effective_component(
  expansion::FloquetExpansion{G,P,E,C,R}, n::Int
) where {G,P,E<:GeneratorComponent,C,R}
  0 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 0:$(expansion.order - 1)"))
  component = reattach(expansion.effective_components[n + 1], expansion.generator.wd, n)::E
  return SQA.simplify(component)::E
end

"""
    micromotion(expansion::FloquetExpansion) -> PeriodicGenerator
    micromotion(expansion::FloquetExpansion, n::Int) -> PeriodicGenerator

Return the retained micromotion generator as harmonics. With `n` in
`1:expansion.order - 1`, return only the order-`n` contribution, including its inverse-frequency
scaling. The micromotion series has no order-0 contribution.

Evaluate the result at symbolic time `t` with `micromotion(expansion)(t)`. Positive completion
leaves the retained micromotion unchanged.

See also [`effective_generator`](@ref), [`effective_component`](@ref), [`VanVleck`](@ref).
"""
function micromotion(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P<:PeriodicGenerator,E,C,R}
  result = zero(expansion.generator)::P
  for (order, kick) in enumerate(expansion.kick_components)
    result = (result + reattach(kick, order))::P
  end
  return result::P
end

function micromotion(
  expansion::FloquetExpansion{G,P,E,C,R}, n::Int
) where {G,P<:PeriodicGenerator,E,C,R}
  1 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 1:$(expansion.order - 1)"))
  return SQA.simplify(reattach(expansion.kick_components[n], n))::P
end
