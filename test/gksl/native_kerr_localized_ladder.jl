include("native_kerr_localized.jl")

# Research-only exact Kerr specialization of the native static homological solve.
#
# For the U(1)-symmetric averaged Kerr generator the charged static map is triangular in the
# two-leg charge level t (native GKLS derivation note, Kerr section):
#
#   y_t = d_t s_t + J s_{t+1}.
#
# The diagonal multiplier d_t is affine in the two occupation numbers, while J is the recycling
# contribution from the static loss channel.  This fixture solves that triangular system with
# exact rational arithmetic on a finite-Fock support.  The finite cutoff only terminates the
# downward charge ladder; every retained coefficient is the same localized/rational coefficient
# that appears in the infinite charged-sector recursion.

const LocalizedCoeff2 = Function

localized_zero_coeff2() = (_, _) -> zero(LocalizedExact)
localized_one_coeff2() = (_, _) -> one(LocalizedExact)

function localized_falling(n, degree)
  degree == 0 && return 1
  return prod(n - offset for offset in 0:(degree - 1))
end

function localized_f_constant_target(qleft, qright)
  left_degree = max(0, -qleft)
  right_degree = max(0, -qright)
  return function (nleft, nright)
    value = localized_falling(nleft, left_degree) * localized_falling(nright, right_degree)
    return LocalizedExact(value, 0)
  end
end

function localized_recycling_value(qleft, qright, above, nleft, nright, κ)
  above === nothing && return zero(LocalizedExact)
  first = (nleft + qleft + 1) * (nright + qright + 1) * above(nleft, nright)
  second = if iszero(nleft) || iszero(nright)
    zero(LocalizedExact)
  else
    nleft * nright * above(nleft - 1, nright - 1)
  end
  return κ * (first - second)
end

function localized_level_solution(qleft, qright, target, above, Δ, χ, κ)
  bohr = localized_kerr_bohr2(qleft, qright, Δ, χ, κ)
  return function (nleft, nright)
    numerator =
      target(nleft, nright) -
      localized_recycling_value(qleft, qright, above, nleft, nright, κ)
    denominator = localized_value(bohr, nleft, nright)
    iszero(denominator) && error("physical Kerr resonance in localized static ladder")
    return numerator / denominator
  end
end

function localized_kerr_static_ladder(
  Q, top, bottom, Δ, χ, κ; top_target=localized_one_coeff2()
)
  bottom <= top || throw(ArgumentError("bottom charge level must not exceed top"))
  solution = Dict{Int,LocalizedCoeff2}()
  above = nothing
  for level in top:-1:bottom
    qleft = level + Q
    qright = level
    target = level == top ? top_target : localized_zero_coeff2()
    current = localized_level_solution(qleft, qright, target, above, Δ, χ, κ)
    solution[level] = current
    above = current
  end
  return solution
end

function localized_exact_shift_superoperator(d, qleft, qright, coefficient)
  result = zeros(LocalizedExact, d^2, d^2)
  for nright in 0:(d - 1), nleft in 0:(d - 1)
    outleft = nleft + qleft
    outright = nright + qright
    0 <= outleft < d || continue
    0 <= outright < d || continue
    input = nleft + 1 + d * nright
    output = outleft + 1 + d * outright
    result[output, input] = coefficient(nleft, nright)
  end
  return result
end

function localized_exact_lindbladian(d, Δ, χ, κ)
  result = zeros(LocalizedExact, d^2, d^2)
  for nright in 0:(d - 1), nleft in 0:(d - 1)
    input = nleft + 1 + d * nright
    result[input, input] +=
      localized_kerr_nojump(nleft, Δ, χ, κ) + conj(localized_kerr_nojump(nright, Δ, χ, κ))
    if nleft >= 1 && nright >= 1
      output = nleft + d * (nright - 1)
      result[output, input] += LocalizedExact(κ * nleft * nright, zero(κ))
    end
  end
  return result
end

function localized_exact_ladder_superoperator(d, Q, solution)
  result = zeros(LocalizedExact, d^2, d^2)
  for level in sort!(collect(keys(solution)); rev=true)
    result += localized_exact_shift_superoperator(d, level + Q, level, solution[level])
  end
  return result
end

function localized_physical_support(d, qleft, qright)
  left0 = max(0, -qleft)
  right0 = max(0, -qright)
  left1 = min(d - 1, d - 1 - qleft)
  right1 = min(d - 1, d - 1 - qright)
  return left0:left1, right0:right1
end

function localized_exact_sector_entries(d, Q)
  entries = Tuple{Int,Int}[]
  for nright in 0:(d - 1), nleft in 0:(d - 1), outright in 0:(d - 1)
    outleft = outright - nright + nleft + Q
    0 <= outleft < d || continue
    push!(entries, (outleft + 1 + d * outright, nleft + 1 + d * nright))
  end
  return entries
end

function localized_exact_sector_solve(L0, target, d, Q)
  entries = localized_exact_sector_entries(d, Q)
  m = length(entries)
  A = zeros(LocalizedExact, m, m)
  for (column, (output, input)) in pairs(entries), (row, (rowout, rowin)) in pairs(entries)
    value = zero(LocalizedExact)
    rowin == input && (value += L0[rowout, output])
    rowout == output && (value -= L0[input, rowin])
    A[row, column] = value
  end
  x = A \ LocalizedExact[target[output, input] for (output, input) in entries]
  result = zeros(LocalizedExact, d^2, d^2)
  for (k, (output, input)) in pairs(entries)
    result[output, input] = x[k]
  end
  # The charge-Q sector must be closed under ad_Lbar for the restricted solve to be the full one.
  L0 * result - result * L0 == target || error("charge sector is not closed under ad_Lbar")
  return result
end

@testset "localized Kerr static ladder solves the full loss homological equation" begin
  Δ = big(1) // big(2)
  χ = big(3) // big(10)
  κ = big(4) // big(5)
  Q = -1
  top = -2
  qtop = (-3, -2)

  # This is the same non-bright charged sector as the first modulated-loss obstruction of the
  # driven Kerr polynomial no-go.  A constant F-representation target carries the falling-factorial leg factors in
  # the faithful c-representation used for the full superoperator commutator.
  topbohr = localized_kerr_bohr2(qtop..., Δ, χ, κ)
  top_target = localized_f_constant_target(qtop...)
  @test topbohr.constant == LocalizedExact(2 // 1, -1 // 1)
  @test topbohr.left == LocalizedExact(0 // 1, 9 // 5)
  @test topbohr.right == LocalizedExact(0 // 1, -6 // 5)

  saved = Dict{Int,Dict{Int,LocalizedCoeff2}}()
  for d in (6, 8, 10)
    # level bottom=-(d-2) is the last sector with physical support for Q=-1.  The next recycling
    # image has left charge -d and therefore vanishes identically in the d-level truncation.
    bottom = -(d - 2)
    solution = localized_kerr_static_ladder(Q, top, bottom, Δ, χ, κ; top_target)
    saved[d] = solution

    for level in top:-1:bottom
      qleft = level + Q
      qright = level
      bohr = localized_kerr_bohr2(qleft, qright, Δ, χ, κ)
      lefts, rights = localized_physical_support(d, qleft, qright)
      @test all(
        !iszero(localized_value(bohr, nleft, nright)) for nleft in lefts, nright in rights
      )
    end

    L0 = localized_exact_lindbladian(d, Δ, χ, κ)
    S = localized_exact_ladder_superoperator(d, Q, solution)
    target = localized_exact_shift_superoperator(d, qtop..., top_target)

    # This is the actual static homological equation [Lbar,S]=Y, including the recycling part
    # of the Lindbladian rather than only the diagonal no-jump Bohr multiplier.
    @test L0 * S - S * L0 == target
  end

  # Increasing the Fock cutoff only extends the downward recycling ladder.  Every superoperator
  # coefficient already present at the smaller cutoff is exactly cutoff independent.
  for (small, large) in ((6, 8), (8, 10)), level in top:-1:(-(small - 2))
    qleft = level + Q
    qright = level
    lefts, rights = localized_physical_support(small, qleft, qright)
    @test !isempty(lefts) && !isempty(rights)
    for nleft in lefts, nright in rights
      @test saved[large][level](nleft, nright) == saved[small][level](nleft, nright)
    end
  end

  # Independent certificate: solve [Lbar,S]=Y by exact elimination over the whole Q=-1 sector
  # without the ladder ansatz.  Every level of this sector has Re d_t = -κ(2t+Q)/2 != 0, so the
  # sector map is injective and the unique solution must be the ladder, including zeros above top.
  let d = 6
    L0 = localized_exact_lindbladian(d, Δ, χ, κ)
    target = localized_exact_shift_superoperator(d, qtop..., top_target)
    independent = localized_exact_sector_solve(L0, target, d, Q)
    @test independent == localized_exact_ladder_superoperator(d, Q, saved[d])
  end

  # The top coefficient is the exact two-leg Schur multiplier from #376 applied to the physical
  # falling-factorial target of this charged sector.
  first = saved[10][top]
  for nleft in 3:8, nright in 2:8
    @test localized_value(topbohr, nleft, nright) * first(nleft, nright) ==
      top_target(nleft, nright)
  end
end
