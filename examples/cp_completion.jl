# # [CP-completion examples](@id cp-completion-examples)
#
# A finite-order Floquet-Lindblad expansion can leave the GKLS cone even when the microscopic
# dynamics is completely positive [Ikeda2021, Schnell2021](@cite). These examples use the same
# analytical models as the validation suite. The driven qubit also compares the algebraic Gram
# construction with the spectral/HCM continuation of Haddadfarshi, Cui, and Mintert
# [Haddadfarshi2015](@cite).

using FloquetExpansions
using Symbolics: @variables

# ## Driven qubit: Gram versus HCM
#
# Consider
#
# ```math
# H(t)=\frac{1}{2}\sigma_z+E\cos(\omega t)\sigma_x,
# ```
#
# with amplitude damping ``\gamma\mathcal D[\sigma_-]``.

pauli = PauliSpace(:cp_doc_qubit)
σx = Pauli(pauli, :sigma, 1)
σy = Pauli(pauli, :sigma, 2)
σz = Pauli(pauli, :sigma, 3)
σminus = (1 // 2) * (σx - im * σy)
σplus = (1 // 2) * (σx + im * σy)

@variables ω::Real t::Real E::Real γ::Real
H = (1 // 2) * σz + E * cos(ω * t) * σx
vv = floquet_expansion(H, ω, t, VanVleck(), 3; channels=(jump(σminus, γ),))

# `order = 3` retains the average and the first two inverse-frequency corrections. With
# ``z=E/\omega``,
#
# ```math
# H_{\mathrm{eff}}^{[2]}=\frac{1}{2}(1-z^2)\sigma_z.
# ```
#
# ### Cartesian Gram completion
#
# `Gram()` does not require a spectral basis.

cartesian = DissipativeFrame(σx, σy, σz)
d_raw = kossakowski(vv, cartesian)
cp_cartesian = positive_completion(vv, Gram(), cartesian)

# The retained Cartesian Kossakowski matrix is
#
# ```math
# d^{[2]}_{xyz}=\frac{\gamma}{4}
# \begin{pmatrix}
# 1 & i(1-z^2) & 0\\
# -i(1-z^2) & 1-2z^2 & 0\\
# 0 & 0 & 2z^2
# \end{pmatrix}.
# ```
#
# ### Adapted frame: Gram and spectral/HCM

adapted = DissipativeFrame(σminus, σplus, σz)
gram = positive_completion(vv, Gram(), adapted)
spectral = positive_completion(vv, Spectral(), adapted)

# The retained target in this frame is
#
# ```math
# d^{[2]}_{\pm z}=\gamma
# \begin{pmatrix}
# 1-z^2 & z^2/2 & 0\\
# z^2/2 & 0 & 0\\
# 0 & 0 & z^2/2
# \end{pmatrix}.
# ```
#
# The coefficient-matched Gram continuation is
#
# ```math
# \widetilde d_{\rm Gram}=\gamma
# \begin{pmatrix}
# 1-z^2+z^4/4 & z^2/2-z^4/4 & 0\\
# z^2/2-z^4/4 & z^4/4 & 0\\
# 0 & 0 & z^2/2
# \end{pmatrix}.
# ```
#
# For the spectral/HCM continuation define
#
# ```math
# A(z)=1-z^2+\frac{z^4}{4},\qquad N(z)=1+\frac{z^4}{4}.
# ```
#
# Then
#
# ```math
# \widetilde d_{\rm HCM}=\frac{\gamma}{N(z)}
# \begin{pmatrix}
# A(z) & \frac{z^2}{2}A(z) & 0\\
# \frac{z^2}{2}A(z) & \frac{z^4}{4}A(z) & 0\\
# 0 & 0 & \frac{z^2}{2}N(z)
# \end{pmatrix}.
# ```
#
# Both are positive and reproduce the same retained coefficients through ``z^2``; their
# difference begins only beyond the retained order.

for n in 0:2
  @assert effective_component(gram, n) == effective_component(vv, n)
  @assert effective_component(spectral, n) == effective_component(vv, n)
end
@assert micromotion(gram) == micromotion(vv)
@assert micromotion(spectral) == micromotion(vv)

# Method-specific data remains behind the common factorization accessor.

gram_data = factorization(gram)
spectral_data = factorization(spectral)
gram_data.onsets
spectral_data.onsets
spectral_data.puiseux

# ## Modulated one- and two-photon loss
#
# Consider two periodically modulated reservoirs,
#
# ```math
# \kappa_1[1+r_1\cos(\omega t)]\mathcal D[a]
# +\kappa_2[1+r_2\cos(\omega t-\pi/2)]\mathcal D[a^2].
# ```
#
# The quarter-cycle phase shift produces
#
# ```math
# [\mathcal D[a],\mathcal D[a^2]]
# =2\mathcal D[a^2]-\mathcal C[a,a^\dagger a^2],
# ```
#
# so the Floquet expansion generates ``a^\dagger a^2`` as a new dark dissipative direction.

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
cp_bosonic = positive_completion(vv_bosonic, Gram(), frame)

# With
#
# ```math
# c=\frac{\kappa_1\kappa_2r_1r_2}{2\omega},
# ```
#
# the retained form is
#
# ```math
# d^{[1]}=\begin{pmatrix}
# \kappa_1 & 0 & c\\
# 0 & \kappa_2-2c & 0\\
# c & 0 & 0
# \end{pmatrix}.
# ```
#
# On the regular stratum ``\kappa_1\ne0`` and ``\kappa_2\ne0``, Gram/Feshbach completion adds
#
# ```math
# \widetilde d-d^{[1]}
# =\begin{pmatrix}
# 0 & 0 & 0\\
# 0 & c^2/\kappa_2 & 0\\
# 0 & 0 & c^2/\kappa_1
# \end{pmatrix}.
# ```
#
# No Fock-space cutoff or Liouville-space vectorization enters this calculation. At a zero of
# either leading rate the dissipative rank changes, so the lower-rank parameter stratum must be
# resolved separately.

channels(cp_bosonic)
factorization(cp_bosonic)

# ## Kerr resonator and symbolic rate conditions

fock_kerr = FockSpace(:cp_doc_kerr)
a_kerr = Destroy(fock_kerr, :a)
J = a_kerr' * a_kerr^2
frame_kerr = DissipativeFrame(J)

@variables ωk::Real tk::Real K::Real γk::Real
H_kerr = K * a_kerr'^2 * a_kerr^2
vv_kerr = floquet_expansion(H_kerr, ωk, tk, VanVleck(), 1; channels=(jump(J, γk),))
cp_kerr = positive_completion(vv_kerr, Gram(), frame_kerr)

positivity_conditions(cp_kerr)
regularity_conditions(cp_kerr)

# The physical rate condition ``\gamma\ge0`` and algebraic nonzero pivot conditions have different
# meanings. The former restricts the physical parameter domain; the latter identifies the
# fixed-rank symbolic stratum used by a particular factorization.
#
# See [Positive completion](@ref positive-completion-manual) for the user API and
# [CP-preserving completion](@ref cp-preserving-completion-theory) for the theory.
#
# ## References
#
# ```@bibliography
# Canonical = false
# ```
