## Reading a symbolic time dependence into harmonics, and writing it back out.

_unwrap(x) = Symbolics.unwrap(x)

_isexpim(x) = Symbolics.iscall(x) && Symbolics.operation(x) === SQA.expim

# `iszero` on a `BasicSymbolic` builds the symbolic equation `0 == 0` rather than returning a
# `Bool`, so it is unusable in a condition. Structural comparison against the literal is.
_issymzero(x) = isequal(Symbolics.value(x), 0)

_dependson(x, v) = any(isequal(Symbolics.value(v)), Symbolics.get_variables(_unwrap(x)))

# Used only to refuse a silent misclassification: a phase buried in a structure `_phase_of` does
# not decompose must raise, never be read as harmonic 0.
function _hasexpim(x)
    y = _unwrap(x)
    _isexpim(y) && return true
    Symbolics.iscall(y) || return false
    return any(_hasexpim, Symbolics.arguments(y))
end

function _addparts(x)
    y = _unwrap(x)
    return if Symbolics.iscall(y) && Symbolics.operation(y) === (+)
        collect(Symbolics.arguments(y))
    else
        [y]
    end
end

# `expim(arg)` is `exp(+i*arg)` with `arg` REAL, while the package convention is
# `exp(-i*m*w*t)`. Split `arg` into `c*w*t + offset` and return `(m, offset) = (-c, offset)`;
# the offset is a constant phase, which belongs with the coefficient rather than the index.
function _harmonic_index(arg, w, t)
    offset = Symbolics.substitute(arg, Dict(w => 0, t => 0))
    c = Symbolics.value(Symbolics.substitute(arg, Dict(w => 1, t => 1)) - offset)

    # Guards `w*t^2` and friends: only a phase linear in `w*t` is periodic at all.
    residual = Symbolics.simplify(arg - (c * w * t + offset))
    _issymzero(residual) || throw(
        ArgumentError(
            "phase $(arg) is not of the form c*$(w)*$(t) + constant, so it is not a whole " *
                "harmonic of the drive",
        ),
    )

    c isa Number || throw(
        ArgumentError("harmonic index $(c) is not a number; it may depend on a symbol other than $(w) and $(t)"),
    )
    m = round(Int, real(c))
    (isapprox(real(c), m; atol = 1.0e-12) && abs(imag(c)) < 1.0e-12) ||
        throw(ArgumentError("harmonic index $(c) is not an integer; the drive is not $(2)π/$(w)-periodic"))
    return -m, offset
end

# `(m, rest, offset)` with `part == rest * expim(offset) * expim(-m*w*t)`. The constant phase is
# returned separately rather than folded into `rest`: `expim` yields an SQA `Coeff`, and neither
# `BasicSymbolic * Coeff` nor `BasicSymbolic * Complex{Num}` exists, so it has to meet the QAdd.
function _phase_of(part, w, t)
    p = _unwrap(part)
    zero_offset = _unwrap(0 * w)

    if _isexpim(p)
        m, offset = _harmonic_index(Symbolics.arguments(p)[1], w, t)
        return m, 1, offset
    end

    if Symbolics.iscall(p)
        op = Symbolics.operation(p)

        if op === (*)
            factors = collect(Symbolics.arguments(p))
            phases = findall(_isexpim, factors)
            length(phases) > 1 &&
                throw(ArgumentError("more than one phase factor in $(part); call `simplify` first"))
            if !isempty(phases)
                m, offset = _harmonic_index(Symbolics.arguments(factors[phases[1]])[1], w, t)
                others = setdiff(eachindex(factors), phases)
                rest = isempty(others) ? 1 : prod(factors[others])
                return m, rest, offset
            end

        elseif op === (/)
            # Reattaching wd puts the phase over a denominator, `exp(-im*t*w)*g / w`. Before this
            # case existed the division fell through and was read as harmonic 0, silently.
            num, den = Symbolics.arguments(p)
            (_hasexpim(den) || _dependson(den, t)) && throw(
                ArgumentError("time dependence in a denominator, in $(part); this is not periodic"),
            )
            m, rest, offset = _phase_of(num, w, t)
            return m, rest / den, offset
        end
    end

    # Anything still carrying a phase here would be misread as the DC harmonic.
    _hasexpim(p) && throw(
        ArgumentError(
            "cannot read a harmonic index from $(part): it carries a phase in a form this " *
                "parser does not decompose. Call `simplify` on the input first.",
        ),
    )
    return 0, p, zero_offset
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

Inverse of calling the result: `harmonics(H, w, t)(w, t)` reproduces `H`.

# Examples

```jldoctest
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
    out = Dict{Int, SQA.QAdd}()
    for (term, coeff) in SQA.exponential_form(H)
        mono = isempty(term.ops) ? one(SQA.QAdd) : prod(term.ops)
        z = SQA.to_num(coeff)
        # A QAdd maps monomial -> coefficient, so ONE coefficient can be a sum over several
        # harmonics. Classifying per term rather than per additive part loses those.
        for (chunk, unit) in ((real(z), 1), (imag(z), im))
            _issymzero(chunk) && continue
            for part in _addparts(chunk)
                m, rest, offset = _phase_of(part, w, t)
                contribution = (unit * rest) * mono
                _issymzero(offset) || (contribution = contribution * SQA.expim(offset))
                out[m] = haskey(out, m) ? out[m] + contribution : contribution
            end
        end
    end
    return PeriodicOperator(out)
end

harmonics(H::SQA.QSym, w, t) = harmonics(_qadd(H), w, t)

"""
    (X::PeriodicOperator)(w, t) -> QAdd

Rebuild the time-dependent operator ``\\sum_l X_l e^{-i l w t}`` from its harmonics, the
inverse of [`harmonics`](@ref).
"""
function (X::PeriodicOperator)(w, t)
    return sum(SQA.expim(-l * w * t) * Xl for (l, Xl) in X.components; init = zero(SQA.QAdd))
end
