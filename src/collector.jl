# `iszero` on a `BasicSymbolic` builds the symbolic equation `0 == 0` rather than returning a
# `Bool`, so use structural comparison instead.
issymzero(x) = isequal(Symbolics.value(x), 0)

# `expim(arg)` is `exp(+i*arg)` with arg REAL; the package convention is `exp(-i*m*w*t)`.
# Split arg into `c*w*t + offset` and return `(-c, offset)`.
function harmonic_index(arg, w, t)
  offset = Symbolics.substitute(arg, Dict(t => 0))
  time_part = Symbolics.simplify(arg - offset)
  c = Symbolics.value(Symbolics.substitute(time_part, Dict(w => 1, t => 1)))

  # Guards `w*t^2` and friends: only a phase linear in `w*t` is periodic at all.
  residual = Symbolics.simplify(arg - (c * w * t + offset))
  issymzero(residual) || throw(
    ArgumentError(
      "phase $(arg) is not of the form c*$(w)*$(t) + constant, so it is not a whole " *
      "harmonic of the drive",
    ),
  )

  c isa Number || throw(
    ArgumentError(
      "harmonic index $(c) is not a number; it may depend on a symbol other than $(w) and $(t)",
    ),
  )
  m = round(Int, real(c))
  (isapprox(real(c), m; atol=1.0e-12) && abs(imag(c)) < 1.0e-12) || throw(
    ArgumentError(
      "harmonic index $(c) is not an integer; the drive is not $(2)π/$(w)-periodic"
    ),
  )
  return -m, offset
end

"""
    harmonics(H::QAdd, w, t) -> PeriodicGenerator

Split a time-dependent operator into its Fourier harmonics, in the convention

```math
H(t) = \\sum_m H_m \\, e^{-i m w t}
```

`w` is the drive frequency and `t` the time variable, both symbolic. Trigonometric time
dependence is normalized to phases first, so `cos`, `sin` and `expim` are all accepted, as is
a constant phase offset such as `cos(w*t + φ)`.

Inverse of calling the result: `harmonics(H, w, t)(t)` reproduces `H`.

# Examples

```jldoctest
julia> using LinearAlgebra: ishermitian

julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real;

julia> H = harmonics(cos(w * t) * (a + a'), w, t)
PeriodicGenerator with harmonics -1:1
  l = -1  =>  1//2 * a + 1//2 * a'
  l = 1  =>  1//2 * a + 1//2 * a'

julia> ishermitian(H)
true
```

See also [`PeriodicGenerator`](@ref).
"""
function harmonics(H::SQA.QAdd, w::Symbolics.Num, t::Symbolics.Num)
  out = Dict{Int,SQA.QAdd}()
  for (term, coeff) in SQA.exponential_form(H)
    mono = isempty(term.ops) ? one(SQA.QAdd) : prod(term.ops)
    for phase_term in SQA.phase_terms(coeff)
      m, offset = harmonic_index(phase_term.phase, w, t)
      contribution = phase_term.amplitude * mono
      issymzero(offset) || (contribution = contribution * SQA.expim(offset))
      out[m] = haskey(out, m) ? out[m] + contribution : contribution
    end
  end
  return PeriodicGenerator(out, w, zero(H))
end

harmonics(H::SQA.QSym, w::Symbolics.Num, t::Symbolics.Num) = harmonics(qadd(H), w, t)

"""
    harmonics(L::Liouvillian, w, t) -> PeriodicGenerator{Liouvillian}

Lower a symbolic time-dependent Liouvillian into the common periodic-generator
representation. The Fourier decomposition is applied independently to the left and right
operator factors of every action and to its scalar coefficient, so periodic dependence in a
Hamiltonian, collapse operator, jump operator, rate, or any combination is supported.

Construct the time-dependent map with [`Liouvillian`](@ref), then pass it here before calling
[`floquet_expansion`](@ref). The result is the same native `PeriodicGenerator{Liouvillian}` as a
manually assembled periodic Liouvillian.
"""
function harmonics(L::Liouvillian, w::Symbolics.Num, t::Symbolics.Num)
  out = Dict{Int,Liouvillian}()
  for ((left, right), coefficient) in term_pairs(L)
    left_harmonics = harmonics(left, w, t)
    right_harmonics = harmonics(right, w, t)

    for (left_harmonic, left_component) in component_pairs(left_harmonics),
      (right_harmonic, right_component) in component_pairs(right_harmonics),
      phase_term in SQA.phase_terms(coefficient)

      phase_coefficient = phase_term.amplitude
      m, offset = harmonic_index(phase_term.phase, w, t)
      issymzero(offset) || (phase_coefficient *= SQA.expim(offset))
      harmonic = left_harmonic + right_harmonic + m

      haskey(out, harmonic) || (out[harmonic] = zero(L))
      add_term!(out[harmonic], left_component, right_component, phase_coefficient)
    end
  end
  return PeriodicGenerator(out, w, zero(L))
end
