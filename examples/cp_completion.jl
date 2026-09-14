# # [CP-preserving completion](@id cp-completion-examples)
#
# A finite-order Floquet expansion of a periodically driven Lindblad generator can leave the
# completely-positive cone. Here we show how to construct a positive completion that preserves the retained Floquet order.
#
# The examples follow the driven-qubit and rank-deficient Gram constructions described in
# [Ikeda2021, Schnell2021](@cite).

using FloquetExpansions

# We start with a qubit whose Hamiltonian is driven along ``x`` and whose microscopic decay is
# the jump ``\sigma_-``.

pauli = PauliSpace(:cp_doc_qubit)
σx = Pauli(pauli, :σ, 1)
σy = Pauli(pauli, :σ, 2)
σz = Pauli(pauli, :σ, 3)
σminus = (1 // 2) * (σx - im * σy)
σplus = (1 // 2) * (σx + im * σy)

@variables ω::Real t::Real E::Real γ::Real
H = (1 // 2) * σz + E * cos(ω * t) * σx

# The Lindblad generator combines the coherent action of ``H`` with the dissipator of
# ``\sigma_-`` at rate ``\gamma``.

L = liouvillian(H; channels=(jump(σminus, γ),))

# We will apply a the floquet expansion to the Lindblad generator using the Van Vleck gauge up to third order:

vv = floquet_expansion(L, ω, t, VanVleck(), 3)

# The effective Hamiltonian contain the RWA term and the Bloch-Siegert shift

hamiltonian(vv)

# The Kossakowski matrix is a representation of the dissipative part of the Lindblad generator in a specific operator basis. We can display the resulting Kossakowski matrix in the Cartesian Pauli basis using the `DissipativeFrame` constructor.

cartesian = DissipativeFrame(σx, σy, σz)
d_raw = kossakowski(vv, cartesian)

# This matrix is not positive semidefinite, hence the truncated generator is not completely
# positive. Indeed its determinant is negative.

using LinearAlgebra
det(d_raw)

# This is detrimental for the physical interpretation of the truncated generator, as it implies that the evolution it generates can lead to negative probabilities. The phenomenon is a consequence of the truncation of the Floquet expansion, which can introduce unphysical artifacts. To remedy this, we can apply a positive completion to the truncated Kossakowski matrix.

# ## Spectral completion

# The spectral/HCM construction of Haddadfarshi, Cui, and Mintert
# [Haddadfarshi2015](@cite) gives the natural first answer to this problem. It follows the
# perturbative decay-rate branches, completes each retained rate through a finite square, and
# reconstructs a positive finite Kossakowski form.
#
# This route needs a suitable *spectral frame*: its leading Kossakowski component must already be
# diagonal, with any degenerate leading sector resolved before branch recursion begins. The
# Cartesian frame above made the defect visible, but does not meet that requirement. For this
# qubit, the jump-adapted frame does.

spectral_frame = DissipativeFrame(σminus, σplus, σz)
d0_spectral = kossakowski_component(vv, spectral_frame, 0)

# The displayed leading matrix has one active decay direction and two dark directions. In this
# frame the spectral construction can follow the rate branches without asking the package to
# discover a symbolic eigenbasis.

spectral = positive_completion(vv, Spectral(), spectral_frame)
d_sprectral = kossakowski(spectral)

# Note that the correction are all of higher order than the retained Floquet expansion. Hence the completed generator is a valid finite-order Floquet expansion of the original Lindblad generator.

# The completed channels make the continuation into a finite GKLS generator. Spectral completion
# naturally presents them as rate-weighted jump channels.

channels(spectral)

# ## Gram completion: preserve the supplied frame
#
# Spectral completion can be annoying as it requires a diagonalization of the leading Kossakowski matrix. Calling `Spectral()` with the Cartesian frame would fail because its
# leading Kossakowski matrix is not diagonal.

# Instead, one can use a different approach by using a `Gram()` factorization. It constructs a graded factor ``B`` such that
# ``\Pi_N(BB^\dagger)=d^{[N]}``. Algebraically, it separates the leading active sector from its
# dark radical, factors the active part, and recursively resolves the Feshbach residual in the
# dark sector. It therefore needs neither a Hilbert-space or Liouville-space matrix nor a
# symbolic eigendecomposition.

gram = positive_completion(vv, Gram(), cartesian)
d_gram = kossakowski(gram)

# As with the spectral result, the higher-order corrections are all of higher order than the retained Floquet expansion.

d_gram - d_raw

# The completed channels are

channels(gram)

# ## Same retained data, different finite completions
#
# Frame choice and algorithm choice are separate. To compare only the algorithms, run Gram in the
# same jump-adapted frame as the spectral construction.

gram_adapted = positive_completion(vv, Gram(), spectral_frame)

# The two finite Kossakowski differ

kossakowski(gram_adapted) - kossakowski(spectral)

# Thus Gram is not a more accurate high-frequency expansion. Its advantage is a broader
# algebraic domain: it avoids spectral preprocessing, works in a supplied physical frame, and
# exposes active/dark onset information when a leading dissipative form is rank deficient.
