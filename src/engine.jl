"""
    FloquetExpansion

Result of [`floquet_expansion`](@ref). Read it with [`effective_generator`](@ref) and
[`micromotion`](@ref) rather than by field access.
"""
struct FloquetExpansion{G<:Gauge,P<:PeriodicGenerator,E}
  generator::P
  kick_components::Vector{P}
  effective_components::Vector{E}
  gauge::G
  order::Int
end

function Base.getproperty(expansion::FloquetExpansion, name::Symbol)
  if name === :kick_derivative_components ||
    name === :dressed_generator ||
    name === :dressed_kick_derivative
    throw(
      ArgumentError(
        "FloquetExpansion field :$(name) is private; use `effective_generator` and `micromotion`",
      ),
    )
  end
  return getfield(expansion, name)
end

function Base.propertynames(::FloquetExpansion{G,P,E}, private::Bool=false) where {G,P,E}
  names = (:generator, :kick_components, :effective_components, :gauge, :order)
  return if private
    (names..., :kick_derivative_components, :dressed_generator, :dressed_kick_derivative)
  else
    names
  end
end

const GeneratorComponent = Union{SQA.QAdd,Liouvillian}

triindex(n::Int, j::Int) = (n * (n + 1)) ÷ 2 + j + 1
weight_generator(j::Int) = im^j * (1 // factorial(j))
weight_kick_derivative(j::Int) = im^j * (1 // factorial(j + 1))

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
) where {P}
  result = zero(generator)
  for j in 0:n
    result = result + weight_generator(j) * dressed_generator[triindex(n, j)]
  end
  for j in 1:n
    result = result - weight_kick_derivative(j) * dressed_kick_derivative[triindex(n, j)]
  end
  return result
end

"""
    floquet_expansion(generator::PeriodicGenerator, gauge::Gauge, order::Int) -> FloquetExpansion
    floquet_expansion(L::Liouvillian, w, t, gauge::Gauge, order::Int) -> FloquetExpansion
    floquet_expansion(H::QField, w, t, gauge::Gauge, order::Int; channels=()) -> FloquetExpansion

Expand a periodically driven generator into a time-independent effective generator and a
periodic micromotion generator, to the given `order` in `1/generator.wd`. The drive frequency is
bound to `generator`.

The second form parses the time dependence first, so `H` may be given as an ordinary
time-dependent Hamiltonian in the drive frequency `w` and time `t`. The `channels` keyword accepts
[`collapse`](@ref) and [`jump`](@ref) channel values and lowers them into a periodic Liouvillian.
The `Liouvillian` form is useful when the native map is constructed separately.

Hamiltonian generators are checked for Hermiticity at ingest, while Liouvillian generators use the
common algebra without a Hamiltonian Hermiticity requirement.

Truncation follows the spec: `order = 1` retains the time average, and the error is
``\\mathcal{O}(w_d^{-\\text{order}})``. Note that the Liouvillian literature counts this shifted
by one.

The series is ASYMPTOTIC, not convergent. Beyond a problem-dependent optimal order, a higher
`order` makes the answer worse rather than better.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real g::Real;

julia> H = harmonics(w * (a' * a) + g * cos(w * t) * (a + a'), w, t);

julia> vv = floquet_expansion(H, VanVleck(), 1)
FloquetExpansion{VanVleck} of order 1

julia> effective_generator(vv)
w * a' * a
```

See also [`effective_generator`](@ref), [`micromotion`](@ref), and [`harmonics`](@ref).
"""
function floquet_expansion(
  generator::P, gauge::G, order::Int
) where {P<:PeriodicGenerator,G<:Gauge}
  order >= 1 || throw(ArgumentError("order must be >= 1"))

  if generator isa PeriodicGenerator{SQA.QAdd}
    LinearAlgebra.ishermitian(generator) || throw(
      ArgumentError(
        "the drive is not Hermitian: it must satisfy H_{-m} = H_m' (eq:fourierH)"
      ),
    )
  end

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

  return FloquetExpansion(generator, K, effective, gauge, order)
end

function floquet_expansion(
  L::Liouvillian, wd::Symbolics.Num, t::Symbolics.Num, gauge::Gauge, order::Int
)
  return floquet_expansion(harmonics(L, wd, t), gauge, order)
end

function floquet_expansion(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  gauge::Gauge,
  order::Int;
  channels::LiouvillianChannelCollection=(),
)
  if isempty(channels)
    return floquet_expansion(harmonics(qadd(H), wd, t), gauge, order)
  end
  L = Liouvillian(H; channels)
  return floquet_expansion(harmonics(L, wd, t), gauge, order)
end

function reattach(component::GeneratorComponent, wd::Symbolics.Num, n::Int)
  return iszero(n) ? component : wd^(-n) * component
end
function reattach(generator::PeriodicGenerator{T}, n::Int) where {T<:GeneratorComponent}
  return iszero(n) ? generator : generator.wd^(-n) * generator
end

"""
    effective_generator(expansion::FloquetExpansion) -> T
    effective_generator(expansion::FloquetExpansion, n::Int) -> T

The time-independent effective generator ``G_\\text{eff}^{[N]} = \\sum_{k<N} G_\\text{eff}^{(k)}``,
or with `n` the order-`n` contribution alone. The drive-frequency scaling is reattached, so the
order-`n` piece carries ``w_d^{-n}``.
"""
function effective_generator(expansion::FloquetExpansion)
  result = zero(expansion.effective_components[1])
  for n in 0:(expansion.order - 1)
    result =
      result + reattach(expansion.effective_components[n + 1], expansion.generator.wd, n)
  end
  return SQA.simplify(result)
end

function effective_generator(expansion::FloquetExpansion, n::Int)
  0 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 0:$(expansion.order - 1)"))
  return SQA.simplify(
    reattach(expansion.effective_components[n + 1], expansion.generator.wd, n)
  )
end

"""
    micromotion(expansion::FloquetExpansion) -> PeriodicGenerator
    micromotion(expansion::FloquetExpansion, n::Int) -> PeriodicGenerator

The periodic micromotion generator ``K^{[N]} = \\sum_{k<N} K^{(k)}``, as harmonics.

With `n` in `1:expansion.order - 1`, the order-`n` contribution alone, frequency-scaled like
[`effective_generator`](@ref). The micromotion series has no order-0 contribution.

Under [`VanVleck`](@ref) this has vanishing time average, which is what makes
[`effective_generator`](@ref) independent of the initial phase of the drive.

Evaluate the returned periodic generator at a symbolic time with `micromotion(expansion)(t)`.
An integer second argument selects an order; use `micromotion(expansion)(t)` for time
evaluation.
"""
function micromotion(expansion::FloquetExpansion)
  result = zero(expansion.generator)
  for (order, kick) in enumerate(expansion.kick_components)
    result = result + reattach(kick, order)
  end
  return result
end

function micromotion(expansion::FloquetExpansion, n::Int)
  1 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 1:$(expansion.order - 1)"))
  return SQA.simplify(reattach(expansion.kick_components[n], n))
end

"""
    effective_hamiltonian(expansion::FloquetExpansion)
    effective_hamiltonian(expansion::FloquetExpansion, n::Int)

Compatibility aliases for [`effective_generator`](@ref) retained for Hamiltonian callers.
"""
effective_hamiltonian(expansion::FloquetExpansion) = effective_generator(expansion)
function effective_hamiltonian(expansion::FloquetExpansion, n::Int)
  return effective_generator(expansion, n)
end

"""
    kick_operator(expansion::FloquetExpansion)
    kick_operator(expansion::FloquetExpansion, n::Int)
    kick_operator(expansion::FloquetExpansion, t::Symbolics.Num)

Compatibility aliases for [`micromotion`](@ref) retained for Hamiltonian callers.
"""
kick_operator(expansion::FloquetExpansion) = micromotion(expansion)
kick_operator(expansion::FloquetExpansion, n::Int) = micromotion(expansion, n)
function kick_operator(expansion::FloquetExpansion, t::Symbolics.Num)
  return micromotion(expansion)(t)
end
