# # [CP-preserving completion](@id cp-completion-examples)
#
# A finite-order Floquet expansion of a periodically driven Lindblad generator can leave the
# completely-positive cone. This tutorial uses Julia's own output to make that phenomenon and
# its completion visible. The prose introduces the physical models; the matrices, channels, and
# factorization data below are all returned by FloquetExpansions.
#
# The examples follow the driven-qubit and rank-deficient Gram constructions described in
# [Ikeda2021, Schnell2021](@cite). They also compare the algebraic Gram completion with the
# spectral/HCM completion of Haddadfarshi, Cui, and Mintert [Haddadfarshi2015](@cite).

using FloquetExpansions

# ## A driven qubit
#
# We start with a qubit whose Hamiltonian is driven along ``x`` and whose microscopic decay is
# the jump ``\sigma_-``. Keeping the symbols unevaluated is useful here: the returned objects show
# the inverse-frequency dependence that the package derives from the Fourier components.

pauli = PauliSpace(:cp_doc_qubit)
σx = Pauli(pauli, :σ, 1)
σy = Pauli(pauli, :σ, 2)
σz = Pauli(pauli, :σ, 3)
σminus = (1 // 2) * (σx - im * σy)
σplus = (1 // 2) * (σx + im * σy)

@variables ω::Real t::Real E::Real γ::Real
H = (1 // 2) * σz + E * cos(ω * t) * σx

# The microscopic model is passed to `floquet_expansion` together with its physical jump channel.
# `order = 3` retains the average and the first two inverse-frequency corrections.

vv = floquet_expansion(H, ω, t, VanVleck(), 3; channels=(jump(σminus, γ),))

# The effective Hamiltonian is available independently of the dissipative representation.

hamiltonian(vv)

# ## Cartesian coordinates: raw versus completed
#
# A `DissipativeFrame` fixes the operator directions in which the Kossakowski form is displayed.
# First use the Cartesian Pauli directions and inspect the raw finite-order result.

cartesian = DissipativeFrame(σx, σy, σz)
d_raw = kossakowski(vv, cartesian)

# `Gram()` constructs a graded collapse-amplitude factor and uses its finite Gram product as the
# positive continuation. The completed Kossakowski matrix therefore contains terms beyond the
# directly retained expansion.

cp_cartesian = positive_completion(vv, Gram(), cartesian)
d_cartesian = kossakowski(cp_cartesian)

# The difference isolates the higher-order information supplied by the completion.

d_cartesian - d_raw

# The completion exposes the corresponding finite collapse channels. Their existence is the
# concrete GKLS representation behind the positive matrix above.

channels(cp_cartesian)

# The factorization accessor gives the graded data used to construct those channels. The onset
# list and the full diagnostic object show that a dark Cartesian direction opens at higher order.

cartesian_factorization = factorization(cp_cartesian)

# The factorization onsets are also available directly.

cartesian_factorization.onsets

# Completion changes only information beyond the retained order. These public-API comparisons
# check every retained Kossakowski component and the micromotion.

[
  kossakowski_component(cp_cartesian, n) == kossakowski_component(vv, cartesian, n)
  for n in 0:2
]

# The micromotion is unchanged as well.

micromotion(cp_cartesian) == micromotion(vv)

# The finite completed generator can be reconstructed from its displayed Hamiltonian and
# collapse channels.

liouvillian(hamiltonian(cp_cartesian); channels=channels(cp_cartesian)) ==
effective_generator(cp_cartesian)

# ## Frame dependence and two completion algorithms
#
# The entries of a Kossakowski matrix depend on the frame. Use the microscopic jump directions as
# an adapted frame, where the leading dissipative form has an active direction and dark directions.

adapted = DissipativeFrame(σminus, σplus, σz)
d_adapted = kossakowski(vv, adapted)

# Both algorithms start from exactly this same retained target.

gram = positive_completion(vv, Gram(), adapted)
spectral = positive_completion(vv, Spectral(), adapted)
kossakowski(gram)

# The spectral continuation is another positive finite matrix in the same frame.

kossakowski(spectral)

# The finite continuations are allowed to differ beyond the retained order, while their retained
# components agree with the raw expansion.

kossakowski(gram) == kossakowski(spectral)

# The Gram continuation agrees with the raw expansion at every retained order.

[
  kossakowski_component(gram, n) == kossakowski_component(vv, adapted, n)
  for n in 0:2
]

# The same coefficient check holds for the spectral continuation.

[
  kossakowski_component(spectral, n) == kossakowski_component(vv, adapted, n)
  for n in 0:2
]

# The two factorization diagnostics make the algorithm-specific information inspectable without
# changing the common `FloquetExpansion` interface.

gram_factorization = factorization(gram)
spectral_factorization = factorization(spectral)
gram_factorization.onsets

# The spectral algorithm records its own onset data.

spectral_factorization.onsets

# Its onset representation also records whether a fractional amplitude is needed.

spectral_factorization.puiseux

# Both completed results preserve the coherent data and reconstruct their effective generators.

micromotion(gram) == micromotion(vv)

# The spectral completion preserves micromotion too.

micromotion(spectral) == micromotion(vv)

# The Gram channels reconstruct its completed generator.

liouvillian(hamiltonian(gram); channels=channels(gram)) == effective_generator(gram)

# The spectral channels reconstruct their completed generator.

liouvillian(hamiltonian(spectral); channels=channels(spectral)) ==
effective_generator(spectral)

# ## Modulated one- and two-photon loss
#
# The next example shows how the package discovers a new dissipative direction generated by the
# Floquet commutator. Two bosonic reservoirs are modulated with a quarter-cycle phase shift.

fock = FockSpace(:cp_doc_bosonic)
a = Destroy(fock, :a)
number_selective = a' * a^2
frame = DissipativeFrame(a, a^2, number_selective)

@variables ωb::Real κ1::Real κ2::Real r1::Real r2::Real
u = κ1 * r1 / 2
v = κ2 * r2 / 2
D1 = dissipator(a)
D2 = dissipator(a^2)

generator = PeriodicGenerator(
  Dict(0 => κ1 * D1 + κ2 * D2, 1 => u * D1 + im * v * D2, -1 => u * D1 - im * v * D2), ωb
)
vv_bosonic = floquet_expansion(generator, VanVleck(), 2)

# The explicit frame lets us inspect the raw form produced by the modulated reservoirs.

kossakowski(vv_bosonic, frame)

# Gram completion adds the higher-order diagonal terms required by the new cross-channel coupling.

cp_bosonic = positive_completion(vv_bosonic, Gram(), frame)
kossakowski(cp_bosonic)

# The completed collapse channels show the new operator direction explicitly.

channels(cp_bosonic)

# The Gram factorization records the same construction in coefficient form.

factorization(cp_bosonic)

# Again, the public accessors verify that the finite channels reconstruct the completed generator.

liouvillian(hamiltonian(cp_bosonic); channels=channels(cp_bosonic)) ==
effective_generator(cp_bosonic)

# ## A rate-modulated qubit
#
# Dissipative Fourier harmonics can mix with coherent harmonics in the same recursion. To expose
# that mechanism directly, modulate the decay rate of a second qubit with a sine phase. The
# opposite imaginary coefficients below are the two nonzero Fourier components of that modulation.

mod_pauli = PauliSpace(:cp_doc_modulated_qubit)
mod_x = Pauli(mod_pauli, :σ, 1)
mod_y = Pauli(mod_pauli, :σ, 2)
mod_z = Pauli(mod_pauli, :σ, 3)
mod_minus = (1 // 2) * (mod_x - im * mod_y)
mod_frame = DissipativeFrame(mod_x, mod_y, mod_z)

@variables ωm::Real tm::Real Δm::Real Em::Real γm::Real am::Real
mod_dissipator = dissipator(mod_minus)
mod_h0 = hamiltonian_action((Δm / 2) * mod_z)
mod_h1 = hamiltonian_action((Em / 2) * mod_x)
mod_rate_harmonic = (γm * am / 2) * mod_dissipator
mod_generator = PeriodicGenerator(
  Dict(
    0 => mod_h0 + γm * mod_dissipator,
    1 => mod_h1 - im * mod_rate_harmonic,
    -1 => mod_h1 + im * mod_rate_harmonic,
  ),
  ωm,
)
mod_vv = floquet_expansion(mod_generator, VanVleck(), 3)

# The first two displayed components show the dissipative targets generated at successive
# inverse-frequency orders.

kossakowski_component(mod_vv, mod_frame, 1)

# The next component shows the order-two dark-sector diffusion.

kossakowski_component(mod_vv, mod_frame, 2)

# The Gram completion succeeds on the symbolic positivity stratum and makes the resulting finite
# channels explicit.

mod_cp = positive_completion(mod_vv, Gram(), mod_frame)
kossakowski(mod_cp)

# The completed collapse channels make the positive continuation constructive rather than merely
# a symbolic claim about a matrix.

channels(mod_cp)

# The factorization object records the active and dark onset data.

factorization(mod_cp)

# The conditions are returned as Julia values, so a caller can inspect the physical assumptions
# and the regularity assumptions used by the symbolic factorization.

positivity_conditions(mod_cp)

# Regularity conditions describe the symbolic stratum used by the factorization.

regularity_conditions(mod_cp)

# As before, the retained coefficients and finite-generator reconstruction are directly testable
# through the public API.

[
  kossakowski_component(mod_cp, n) == kossakowski_component(mod_vv, mod_frame, n)
  for n in 0:2
]

# The finite channels reconstruct the completed effective generator.

liouvillian(hamiltonian(mod_cp); channels=channels(mod_cp)) ==
effective_generator(mod_cp)

# The package therefore exposes the full workflow as inspectable Julia values: symbolic effective
# generators, frame-dependent Kossakowski coordinates, positive finite completions, collapse
# channels, factorization diagnostics, and the assumptions needed by symbolic elimination.
#
# ## References
#
# ```@bibliography
# Canonical = false
# ```
