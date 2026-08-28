using Test
using FloquetExpansions
import SecondQuantizedAlgebra as SQA
using Symbolics: Symbolics

h = FockSpace(:cavity)
a = Destroy(h, :a)
@variables w::Real g::Real ξ::Real t::Real

# Harmonics -2..2 with different operator content each, so the oracle sums below cannot agree
# by degeneracy.
function drive()
  return PeriodicOperator(
    0 => 1 * (a' * a),
    1 => g * a,
    -1 => conj(g) * a',
    2 => ξ * (a' * a'),
    -2 => conj(ξ) * (a * a),
  )
end

comm = SQA.commutator
vanishes(q::SQA.QAdd) = iszero(SQA.simplify(q))

# Exact comparison is the rule (DESIGN D5), but it is not always reachable: SQA's `Native`
# coefficient tier is ComplexF64 and swallows any float-representable rational, so `1//2` is
# stored as `0.5` and every later product with it is done in floating point. A chain of exact
# rationals stays exact ((1//3)(1//5)(1//7)(1//9)(1//11) -> 1//10395); insert one `1//2` and the
# same chain lands on 4.81e-5. Where a weight chain has passed through such a rational, compare
# numerically instead of pretending the residue is not there.
function maxcoeff(q::SQA.QAdd, subs)
  m = 0.0
  for (_, coeff) in q
    z = SQA.to_num(coeff)
    re = Symbolics.value(Symbolics.substitute(real(z), subs))
    im_ = Symbolics.value(Symbolics.substitute(imag(z), subs))
    m = max(m, abs(complex(Float64(re), Float64(im_))))
  end
  return m
end

const SUBS = Dict(w => 1.0, g => 0.7, ξ => 0.3, t => 0.4)
function agrees(X::PeriodicOperator, Y::PeriodicOperator)
  return all(
    maxcoeff(SQA.simplify(X[l] - Y[l]), SUBS) < 1.0e-12 for l in union(keys(X), keys(Y))
  )
end

@testset "the spec's closed forms, orders 0 to 2" begin
  # The only independent oracle that exists: sec:vv writes these out in full, so they check the
  # recursion itself rather than checking it against another run of itself. All EXACT.
  H = drive()
  Ms = [-2, -1, 1, 2]
  vv = floquet_expansion(H, w, VanVleck(), 3)

  @test vanishes(effective_hamiltonian(vv, 0) - H[0])                          # eq:K1

  K1 = FloquetExpansions._reattach(vv.K[1], vv.wd, 1)                          # eq:K1
  @test sort!(collect(keys(K1))) == Ms
  @test all(vanishes(K1[m] - (im / (m * w)) * H[m]) for m in Ms)

  # eq:Heff1  =  -1/2 sum_{m!=0} [H_m, H_-m] / (m wd)
  @test vanishes(
    effective_hamiltonian(vv, 1) -
    sum((-1 // 2) * (1 / (m * w)) * comm(H[m], H[-m]) for m in Ms),
  )

  # eq:K2  =  -i sum_m [ [H_m,H_0]/(m^2 wd^2) + 1/2 sum_{m'!=0,m} [H_m',H_{m-m'}]/(m m' wd^2) ]
  K2 = FloquetExpansions._reattach(vv.K[2], vv.wd, 2)
  function oracleK2(m)
    acc = comm(H[m], H[0]) * (1 / (m^2 * w^2))
    for mp in Ms
      mp == m && continue
      acc = acc + comm(H[mp], H[m - mp]) * (1 / (2 * m * mp * w^2))
    end
    return -im * acc
  end
  @test !haskey(K2.components, 0)                                              # van Vleck gauge
  @test all(vanishes(K2[m] - oracleK2(m)) for m in -4:4 if m != 0)

  # eq:Heff2  =  sum_m [[H_-m,H_0],H_m]/(2 m^2 wd^2)
  #            + sum_m sum_{m'!=0,m} [[H_-m,H_{m-m'}],H_m']/(3 m m' wd^2)
  oracleHeff2 = zero(SQA.QAdd)
  for m in Ms
    oracleHeff2 = oracleHeff2 + comm(comm(H[-m], H[0]), H[m]) * (1 / (2 * m^2 * w^2))
    for mp in Ms
      mp == m && continue
      oracleHeff2 =
        oracleHeff2 + comm(comm(H[-m], H[m - mp]), H[mp]) * (1 / (3 * m * mp * w^2))
    end
  end
  @test vanishes(effective_hamiltonian(vv, 2) - oracleHeff2)
end

@testset "the peeling recursion equals its definition as composition sums" begin
  # The triangle is an O(n^2) rewrite of exponential-size index sums. Checking it against those
  # sums directly is what validates the rewrite: the k-range `1:n-j+1` and the Deprit weights
  # are exactly what a naive derivation gets wrong, and the spec's closed forms stop at order 2.
  H = drive()
  N = 4
  vv = floquet_expansion(H, w, VanVleck(), N)

  # ordered tuples (k_1..k_j), each >= 1, summing to n
  function compositions(n::Int, j::Int)
    j == 0 && return n == 0 ? [Int[]] : Vector{Int}[]
    out = Vector{Int}[]
    for first in 1:(n - j + 1), rest in compositions(n - first, j - 1)
      push!(out, [first; rest])
    end
    return out
  end

  # dressedH^(n)_[j] := (i^j/j!) sum_{k_1+..+k_j=n} ad_{K^(k_1)}..ad_{K^(k_j)} H_S
  function refH(n, j)
    acc = zero(PeriodicOperator)
    for ks in compositions(n, j)
      term = H
      for k in Iterators.reverse(ks)          # ad_{K^(k_1)} is the OUTERMOST bracket
        term = comm(vv.K[k], term)
      end
      acc = acc + term
    end
    return (im^j * (1 // factorial(j))) * acc
  end

  # dressedKdot^(n)_[j] := (i^j/(j+1)!) sum_{k_0+..+k_j=n+1} ad_{K^(k_1)}..ad_{K^(k_j)} Kdot^(k_0)
  function refKdot(n, j)
    acc = zero(PeriodicOperator)
    for ks in compositions(n + 1, j + 1)
      term = vv.Kdot[ks[1]]
      for k in Iterators.reverse(ks[2:end])
        term = comm(vv.K[k], term)
      end
      acc = acc + term
    end
    return (im^j * (1 // factorial(j + 1))) * acc
  end

  for n in 0:(N - 1), j in 0:n
    node = FloquetExpansions._weightH(j) * vv.dressedH[FloquetExpansions._tri(n, j)]
    @test agrees(node, refH(n, j))
    if j >= 1
      nodeK =
        FloquetExpansions._weightKdot(j) * vv.dressedKdot[FloquetExpansions._tri(n, j)]
      @test agrees(nodeK, refKdot(n, j))
    end
  end
end

@testset "truncation follows the spec, X^[N] = sum_{k<N}" begin
  H = drive()

  # order 1 is the rotating-wave approximation: H_0 alone, and no micromotion at all.
  rwa = floquet_expansion(H, w, VanVleck(), 1)
  @test vanishes(effective_hamiltonian(rwa) - H[0])
  @test iszero(kick_operator(rwa))

  # K^[N] excludes K^(N), so raising the order by one adds exactly one kick order and one Heff
  # order, and leaves the lower ones untouched.
  for N in 2:4
    lo = floquet_expansion(H, w, VanVleck(), N - 1)
    hi = floquet_expansion(H, w, VanVleck(), N)
    @test all(
      vanishes(effective_hamiltonian(hi, n) - effective_hamiltonian(lo, n)) for
      n in 0:(N - 2)
    )
    @test agrees(
      kick_operator(hi) - kick_operator(lo),
      FloquetExpansions._reattach(hi.K[N - 1], hi.wd, N - 1),
    )
  end
end

@testset "the kick round-trips through the collector" begin
  # This is what caught the collector's silent division bug: reattaching wd puts the phase over
  # a denominator, which used to parse as harmonic 0 without complaint.
  vv = floquet_expansion(drive(), w, VanVleck(), 3)
  K, back = kick_operator(vv), harmonics(kick_operator(vv, t), w, t)
  @test sort!(collect(keys(back))) == sort!(collect(keys(K)))
  @test all(vanishes(back[l] - K[l]) for l in keys(K))
end

@testset "a non-Hermitian drive is refused at ingest" begin
  # H_{-m} = H_m' (eq:fourierH). Without this the expansion runs and silently produces a
  # non-Hermitian effective Hamiltonian at order 3.
  bad = PeriodicOperator(1 => a, -1 => a)
  @test_throws ArgumentError floquet_expansion(bad, w, VanVleck(), 2)
  @test_throws ArgumentError floquet_expansion(drive(), w, VanVleck(), 0)
end

@testset "the parsing entry point agrees with the harmonic one" begin
  H = 1 * (a' * a) + g * cos(w * t) * (a + a')
  viaparse = floquet_expansion(H, w, t, VanVleck(), 3)
  viaharm = floquet_expansion(harmonics(H, w, t), w, VanVleck(), 3)
  @test vanishes(effective_hamiltonian(viaparse) - effective_hamiltonian(viaharm))
end
