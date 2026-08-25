using SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra
using SymbolicUtils, Symbolics
using Test

unw = Symbolics.unwrap
isexpim(x) = SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === SQA.expim

addparts(x) = (x = unw(x);
    SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === (+) ?
        collect(SymbolicUtils.arguments(x)) : [x])

"""(m, rest) with part == rest * expim(-m*w*t).  Our convention: H_S(t)= sum_m H_m exp(-i m w t)."""
function phase_of(part, w, t)
    part = unw(part)
    getm(arg) = begin
        c = Symbolics.value(Symbolics.substitute(arg, Dict(w => 1, t => 1)))
        mi = round(Int, real(c))
        (isapprox(real(c), mi; atol = 1e-12) && abs(imag(c)) < 1e-12) ||
            error("non-integer harmonic exponent: $c")
        -mi                                     # expim(c*w*t) == exp(+i c w t) == exp(-i(-c) w t)
    end
    isexpim(part) && return (getm(SymbolicUtils.arguments(part)[1]), 1)
    if SymbolicUtils.iscall(part) && SymbolicUtils.operation(part) === (*)
        fs = collect(SymbolicUtils.arguments(part))
        idx = findall(isexpim, fs)
        length(idx) > 1 && error("more than one phase factor in $part")
        isempty(idx) && return (0, part)
        rest = isempty(setdiff(eachindex(fs), idx)) ? 1 :
               prod(fs[setdiff(eachindex(fs), idx)])
        return (getm(SymbolicUtils.arguments(fs[idx[1]])[1]), rest)
    end
    (0, part)
end

function harmonics(q::SQA.QAdd, w, t)
    out = Dict{Int, SQA.QAdd}()
    push!(o, l, v) = (out[l] = haskey(out, l) ? out[l] + v : v)
    for (term, coeff) in q
        mono = isempty(term.ops) ? one(SQA.QAdd) : prod(term.ops)
        z = SQA.to_num(coeff)
        for (chunk, unit) in ((real(z), 1), (imag(z), im))
            c0 = unw(chunk); (c0 isa Number && iszero(c0)) && continue
            for part in addparts(chunk)
                m, rest = phase_of(part, w, t)
                push!(out, m, (unit * rest) * mono)
            end
        end
    end
    Dict(l => v for (l, v) in out if !iszero(v))
end

# ---------------- test ----------------
h = FockSpace(:c); a = Destroy(h, :a)
@variables w t g ξ

H = (3//4)*(a'*a) +
    (1//3)*(a'*a'*a*a) +
    g*SQA.expim(-1*w*t)*(a'*a*a) +
    conj(g)*SQA.expim(1*w*t)*(a'*a'*a) +
    ξ*SQA.expim(-2*w*t)*(a'*a') +
    2*SQA.expim(-1*w*t)*a +
    5*SQA.expim(2*w*t)*a

println("input H(t) = ", H, "\n")
hs = harmonics(H, w, t)
for m in sort(collect(keys(hs))); println("  m = ", lpad(m,2), " :  ", hs[m]); end

# round trip
println("\nround-trip check:")
rebuilt = sum(SQA.expim(-m*w*t)*v for (m,v) in hs)
diff = SQA.simplify(rebuilt - H)
println("  simplify(rebuilt - H) = ", diff, "   iszero = ", iszero(diff))
@test iszero(diff)
println("\nexactness: coefficient of m=0 is ", hs[0], "  (must show 3//4, not 0.75)")
