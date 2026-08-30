## Reading a symbolic time dependence into harmonics, and writing it back out.

# `iszero` on a `BasicSymbolic` builds the symbolic equation `0 == 0` rather than returning a
# `Bool`, so it is unusable in a condition. Structural comparison against the literal is.
_issymzero(x) = isequal(Symbolics.value(x), 0)

# `expim(arg)` is `exp(+i*arg)` with `arg` REAL, while the package convention is
# `exp(-i*m*w*t)`. Split `arg` into `c*w*t + offset` and return `(m, offset) = (-c, offset)`;
# the offset is a constant phase, which belongs with the coefficient rather than the index.
function _harmonic_index(arg, w, t)
  offset = Symbolics.substitute(arg, Dict(t => 0))
  time_part = Symbolics.simplify(arg - offset)
  c = Symbolics.value(Symbolics.substitute(time_part, Dict(w => 1, t => 1)))

  # Guards `w*t^2` and friends: only a phase linear in `w*t` is periodic at all.
  residual = Symbolics.simplify(arg - (c * w * t + offset))
  _issymzero(residual) || throw(
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
    harmonics(H::QAdd, w, t) -> PeriodicOperator

Split a time-dependent operator into its Fourier harmonics, in the convention

```math
H(t) = \\sum_m H_m \\, e^{-i m w t}
```

`w` is the drive frequency and `t` the time variable, both symbolic. Trigonometric time
dependence is normalized to phases first, so `cos`, `sin` and `expim` are all accepted, as is
a constant phase offset such as `cos(w*t + φ)`.

Throws if a phase is not of the form `c*w*t + constant` with integer `c`, i.e. if the operator
is not periodic at the drive frequency.

Inverse of calling the result: `harmonics(H, w, t)(t)` reproduces `H`.

# Examples

```jldoctest
julia> using LinearAlgebra: ishermitian

julia> h = FockSpace(:cavity); a = Destroy(h, :a);

julia> @variables w::Real t::Real;

julia> H = harmonics(cos(w * t) * (a + a'), w, t)
PeriodicOperator with harmonics -1:1
  l = -1  =>  0.5 * a + 0.5 * a'
  l =  1  =>  0.5 * a + 0.5 * a'

julia> ishermitian(H)
true
```

See also [`PeriodicOperator`](@ref).
"""
function harmonics(H::SQA.QAdd, w, t)
  out = Dict{Int,SQA.QAdd}()
  for (term, coeff) in SQA.exponential_form(H)
    mono = isempty(term.ops) ? one(SQA.QAdd) : prod(term.ops)
    for phase_term in SQA.phase_terms(coeff)
      m, offset = _harmonic_index(phase_term.phase, w, t)
      contribution = phase_term.amplitude * mono
      _issymzero(offset) || (contribution = contribution * SQA.expim(offset))
      out[m] = haskey(out, m) ? out[m] + contribution : contribution
    end
  end
  return PeriodicOperator(out, w)
end

harmonics(H::SQA.QSym, w, t) = harmonics(_qadd(H), w, t)

"""
    (X::PeriodicOperator)(t) -> QAdd

Rebuild the time-dependent operator ``\\sum_l X_l e^{-i l \\omega_d t}`` from its harmonics, the
inverse of [`harmonics`](@ref). The drive frequency is taken from `X.wd`.
"""
function (X::PeriodicOperator)(t)
  return sum(SQA.expim(-l * X.wd * t) * Xl for (l, Xl) in X.components; init=zero(SQA.QAdd))
end
