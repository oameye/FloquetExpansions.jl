"""
    harmonic_balance(H, embedding::CarrierEmbedding) -> QAdd
    harmonic_balance(L::Liouvillian, embedding::CarrierEmbedding) -> Liouvillian

Construct the unreduced quantum-harmonic-balance generator of a physical Hamiltonian `H` or
Liouvillian `L` in the carrier modes of `embedding`:

```math
H_{\\mathrm{QHB}} = λ_H \\, P\\big[I_t(H)\\big] + H_{\\mathrm{fr}}, \\qquad
H_{\\mathrm{fr}} = -\\sum_j ω_j \\left(a_{j+}^\\dagger a_{j+} - a_{j-}^\\dagger a_{j-}\\right).
```

``I_t`` is the carrier embedding, ``λ_H = 2N_h``, and ``P`` keeps the processes whose total
frequency vanishes identically. For a Liouvillian, ``P`` acts on each ``ρ ↦ AρB`` term as a
whole, so a single physical dissipator distributes its rate over the carriers; the frame
term enters as its Hamiltonian action. The result is time independent.

For a Hamiltonian quadratic in `a` and `a†`, together with its linear drives and dissipators
in `a` and `a†`, the carrier dynamics reconstructs the physical dynamics exactly through
[`reconstruct`](@ref), including counter-rotating terms. For a nonlinear Hamiltonian the
projection ``P`` is an approximation.

Time dependence of `H` must be written in unit phases `expim(ν * t)` or trigonometric
functions of `t` with frequencies built from the same symbols as the carriers and exact
coefficients; any other dependence on `t` throws an `ArgumentError`.

# Examples

```jldoctest
julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables ω₀::Real ω::Real t::Real;

julia> embedding = CarrierEmbedding(a, [ω]; t);

julia> a₊, a₋ = only(carrier_modes(embedding));

julia> H = harmonic_balance(ω₀ * a' * a, embedding);

julia> isequal(simplify(H - ((ω₀ - ω) * a₊' * a₊ + (ω₀ + ω) * a₋' * a₋)), zero(H))
true
```

See also [`CarrierEmbedding`](@ref), [`reconstruct`](@ref), [`carrier_modes`](@ref).
"""
function harmonic_balance(H::SQA.QField, embedding::CarrierEmbedding)
  embedded = reconstruct(embedding, H)
  resonant = static_part(embedded, embedding.time)
  return normalization(embedding) * resonant + frame_hamiltonian(embedding)
end

function harmonic_balance(L::Liouvillian, embedding::CarrierEmbedding)
  t = embedding.time
  resonant = zero(L)
  for ((left, right), coefficient) in term_pairs(L)
    left_components = frequency_components(reconstruct(embedding, left), t)
    right_components = frequency_components(reconstruct(embedding, right), t)
    scalar_components = frequency_components(coefficient, t)
    for (left_frequency, left_component) in left_components,
      (right_frequency, right_component) in right_components,
      (scalar_frequency, scalar_component) in scalar_components

      total = Symbolics.expand(left_frequency + right_frequency + scalar_frequency)
      issymzero(total) || continue
      add_term!(resonant, left_component, right_component, scalar_component)
    end
  end
  return scale(normalization(embedding), resonant) +
         hamiltonian_action(frame_hamiltonian(embedding))
end

normalization(embedding::CarrierEmbedding) = 2 * length(embedding.frequencies)

function frame_hamiltonian(embedding::CarrierEmbedding)
  return sum(
    -ω * (plus' * plus - minus' * minus) for
    (ω, (plus, minus)) in zip(embedding.frequencies, embedding.modes)
  )
end

function depends_on_time(x::Symbolics.Num, t::Symbolics.Num)::Bool
  return any(isequal(Symbolics.value(t)), Symbolics.get_variables(x))
end

function depends_on_time(x::SQA.CNum, t::Symbolics.Num)::Bool
  value = SQA.to_num(x)
  return depends_on_time(real(value), t) || depends_on_time(imag(value), t)
end

function phase_frequency(phase::Symbolics.Num, t::Symbolics.Num)
  frequency = Symbolics.expand(Symbolics.derivative(phase, t))
  depends_on_time(frequency, t) &&
    throw(ArgumentError("phase $(phase) is not linear in the time $(t)"))
  validate_exact_frequency(frequency)
  offset = Symbolics.substitute(phase, Dict(t => 0))
  return frequency, offset
end

function phase_components(coefficient::SQA.CNum, t::Symbolics.Num)
  components = Tuple{Symbolics.Num,SQA.CNum}[]
  for phase_term in SQA.phase_terms(coefficient)
    amplitude = phase_term.amplitude
    depends_on_time(amplitude, t) && throw(
      ArgumentError(
        "coefficient $(coefficient) depends on the time $(t) other than through a phase"
      ),
    )
    frequency, offset = phase_frequency(Symbolics.Num(phase_term.phase), t)
    issymzero(offset) || (amplitude *= SQA.expim(offset))
    push!(components, (frequency, amplitude))
  end
  return components
end

function frequency_components(coefficient::SQA.CNum, t::Symbolics.Num)
  return phase_components(coefficient, t)
end

function frequency_components(operator::SQA.QAdd, t::Symbolics.Num)
  grouped = Dict{Symbolics.Num,SQA.QAdd}()
  for (term, coefficient) in SQA.exponential_form(operator)
    monomial = isempty(term.ops) ? one(SQA.QAdd) : prod(term.ops)
    for (frequency, amplitude) in phase_components(coefficient, t)
      contribution = amplitude * monomial
      grouped[frequency] =
        haskey(grouped, frequency) ? grouped[frequency] + contribution : contribution
    end
  end
  return collect(grouped)
end

function static_part(operator::SQA.QAdd, t::Symbolics.Num)
  result = zero(operator)
  for (frequency, component) in frequency_components(operator, t)
    issymzero(frequency) && (result += component)
  end
  return result
end
