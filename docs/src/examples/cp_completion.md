```@meta
CurrentModule = FloquetExpansions
CollapsedDocStrings = true
```

# CP-completion examples

These examples use the public API and the same analytical reference models used by the physical
validation suite. The purpose is not merely to show syntax: each example illustrates a distinct
piece of the CP-completion physics.

- the driven qubit shows direct Cartesian [`Gram`](@ref) completion and the independent
  [`Spectral`](@ref)/HCM continuation in an adapted frame;
- the modulated bosonic reservoir shows a Floquet-generated dark dissipative direction and its
  Feshbach/Gram closure without a Hilbert-space cutoff;
- the Kerr resonator shows how symbolic rate assumptions are exposed through the common API.

## Driven qubit: Gram versus HCM

Consider

```math
H(t)=\frac{1}{2}\sigma_z+E\cos(\omega t)\sigma_x,
```

with amplitude damping ``\gamma\mathcal D[\sigma_-]``. In Julia,

```julia
using FloquetExpansions
using Symbolics: @variables

pauli = PauliSpace(:cp_doc_qubit)
σx = Pauli(pauli, :sigma, 1)
σy = Pauli(pauli, :sigma, 2)
σz = Pauli(pauli, :sigma, 3)
σminus = (1 // 2) * (σx - im * σy)
σplus = (1 // 2) * (σx + im * σy)

@variables ω::Real t::Real E::Real γ::Real
H = (1 // 2) * σz + E * cos(ω * t) * σx
vv = floquet_expansion(H, ω, t, VanVleck(), 3; channels=(jump(σminus, γ),))
```

`order = 3` retains the average, first inverse-frequency term, and second inverse-frequency term.
With ``z=E/\omega``, the coherent sector is

```math
H_{\mathrm{eff}}^{[2]}
=\frac{1}{2}(1-z^2)\sigma_z,
```

so the second-order correction is ``-E^2\sigma_z/(2\omega^2)``.

### Cartesian Gram completion needs no spectral basis

Use the Cartesian dissipative frame directly:

```julia
cartesian = DissipativeFrame(σx, σy, σz)
d_raw = kossakowski(vv, cartesian)
cp_cartesian = positive_completion(vv, Gram(), cartesian)
```

The retained Kossakowski matrix is

```math
d^{[2]}_{xyz}
=\frac{\gamma}{4}
\begin{pmatrix}
1 & i(1-z^2) & 0\\
-i(1-z^2) & 1-2z^2 & 0\\
0 & 0 & 2z^2
\end{pmatrix}.
```

The leading Cartesian block is not being converted to an eigenbasis before completion. `Gram()`
works with the Hermitian form as supplied, using congruence elimination, graded ``LDL^\dagger``, and
the active/dark construction. This is the characteristic non-spectral use case of the algebraic
method.

### Adapted frame for a direct Gram--HCM comparison

For an explicit decay-rate comparison it is convenient to use

```julia
adapted = DissipativeFrame(σminus, σplus, σz)
gram = positive_completion(vv, Gram(), adapted)
spectral = positive_completion(vv, Spectral(), adapted)
```

This is the analytical ``(u,v,e_z)`` frame used in the derivation, with ``u`` and ``v`` rescaled to
``\sigma_-`` and ``\sigma_+`` so that the public frame contains no symbolic ``\sqrt 2`` factors.
The retained target is

```math
d^{[2]}_{\pm z}
=\gamma
\begin{pmatrix}
1-z^2 & z^2/2 & 0\\
z^2/2 & 0 & 0\\
0 & 0 & z^2/2
\end{pmatrix}.
```

The coefficient-matched Gram/root continuation is

```math
\widetilde d_{\rm Gram}
=\gamma
\begin{pmatrix}
1-z^2+z^4/4 & z^2/2-z^4/4 & 0\\
z^2/2-z^4/4 & z^4/4 & 0\\
0 & 0 & z^2/2
\end{pmatrix}.
```

The spectral/HCM continuation instead follows the perturbative decay-rate branch. Define

```math
A(z)=1-z^2+\frac{z^4}{4},
\qquad
N(z)=1+\frac{z^4}{4}.
```

Then

```math
\widetilde d_{\rm HCM}
=\frac{\gamma}{N(z)}
\begin{pmatrix}
A(z) & \frac{z^2}{2}A(z) & 0\\
\frac{z^2}{2}A(z) & \frac{z^4}{4}A(z) & 0\\
0 & 0 & \frac{z^2}{2}N(z)
\end{pmatrix}.
```

Both finite matrices are positive continuations and reproduce the same retained Kossakowski
coefficients through ``z^2``. They are nevertheless different finite functions: their difference
starts at ``\mathcal O(z^4)``, i.e. beyond the retained order. This is completion nonuniqueness in a
concrete physical model, not a disagreement about the perturbative Floquet data.

The retained dynamics and micromotion can be checked through the common API:

```julia
for n in 0:2
  effective_component(gram, n) == effective_component(vv, n)
  effective_component(spectral, n) == effective_component(vv, n)
end

micromotion(gram) == micromotion(vv)
micromotion(spectral) == micromotion(vv)

liouvillian(hamiltonian(gram); channels=channels(gram)) == effective_generator(gram)
liouvillian(hamiltonian(spectral); channels=channels(spectral)) == effective_generator(spectral)
```

Algorithm-specific information remains behind [`factorization`](@ref):

```julia
gram_data = factorization(gram)
spectral_data = factorization(spectral)

gram_data.onsets
spectral_data.onsets
spectral_data.puiseux
```

## Modulated one- and two-photon loss

Now consider a bosonic mode with two periodically modulated reservoirs,

```math
\kappa_1[1+r_1\cos(\omega t)]\,\mathcal D[a]
+\kappa_2[1+r_2\cos(\omega t-\pi/2)]\,\mathcal D[a^2].
```

The quarter-cycle phase shift makes the first Floquet correction proportional to the dissipator
commutator. In symbolic second-quantized form,

```math
[\mathcal D[a],\mathcal D[a^2]]
=2\mathcal D[a^2]-\mathcal C[a,a^\dagger a^2],
```

so the Floquet expansion generates the new dissipative direction ``a^\dagger a^2`` even though it
was absent from the microscopic channel list.

Construct the harmonics directly:

```julia
fock = FockSpace(:cp_doc_bosonic)
a = Destroy(fock, :a)
number_selective = a' * a^2
frame = DissipativeFrame(a, a^2, number_selective)

@variables ω::Real κ1::Real κ2::Real r1::Real r2::Real
u = κ1 * r1 / 2
v = κ2 * r2 / 2
D1 = dissipator(a)
D2 = dissipator(a^2)

generator = PeriodicGenerator(
  Dict(
    0 => κ1 * D1 + κ2 * D2,
    1 => u * D1 + im * v * D2,
    -1 => u * D1 - im * v * D2,
  ),
  ω,
)

vv = floquet_expansion(generator, VanVleck(), 2)
cp = positive_completion(vv, Gram(), frame)
```

Define

```math
c=\frac{\kappa_1\kappa_2 r_1r_2}{2\omega}.
```

Through the retained first inverse-frequency order,

```math
d^{[1]}
=
\begin{pmatrix}
\kappa_1 & 0 & c\\
0 & \kappa_2-2c & 0\\
c & 0 & 0
\end{pmatrix}.
```

The ``a^\dagger a^2`` direction is dark at leading order but couples to the active one-photon
sector at order ``1/\omega``. Positivity therefore requires a dark-sector contribution at the next
order. On the regular stratum ``\kappa_1\ne0`` and ``\kappa_2\ne0``, the Gram continuation is

```math
\widetilde d
=
\begin{pmatrix}
\kappa_1 & 0 & c\\
0 & \kappa_2-2c+c^2/\kappa_2 & 0\\
c & 0 & c^2/\kappa_1
\end{pmatrix},
```

so the explicitly added higher-order closure is

```math
\widetilde d-d^{[1]}
=
\begin{pmatrix}
0 & 0 & 0\\
0 & c^2/\kappa_2 & 0\\
0 & 0 & c^2/\kappa_1
\end{pmatrix}.
```

At ``\kappa_1=0`` or ``\kappa_2=0`` the leading dissipative rank changes and the corresponding
lower-rank stratum must be treated separately rather than by continuing the displayed rational
factorization through the singular point.

This is the Feshbach/Schur closure of the retained active--dark coupling. The construction is done
entirely in symbolic operator algebra: no Fock-space cutoff, finite Hilbert-space matrix, or
Liouville vectorization is introduced.

That algebraic statement should not be confused with a functional-analytic theorem on the full
bosonic Fock space. For unbounded jump operators, rigorous semigroup generation and domain questions
require a separate analysis beyond the symbolic finite-frame completion.

The completed jump representation is available directly:

```julia
channels(cp)
factorization(cp)
```

and reconstructs the same finite completed generator together with `hamiltonian(cp)`.

## Kerr resonator and symbolic rate conditions

A static number-selective loss channel already illustrates the distinction between physical rate
conditions and algebraic completion diagnostics:

```julia
fock = FockSpace(:cp_doc_kerr)
a = Destroy(fock, :a)
J = a' * a^2
frame = DissipativeFrame(J)

@variables ω::Real t::Real K::Real γ::Real
H = K * a'^2 * a^2
vv = floquet_expansion(H, ω, t, VanVleck(), 1; channels=(jump(J, γ),))

cp = positive_completion(vv, Gram(), frame)
positivity_conditions(cp)
regularity_conditions(cp)
```

The physical channel requires ``\gamma\ge0``; that assumption belongs to
[`positivity_conditions`](@ref). [`regularity_conditions`](@ref) has a different role: it records
nonzero pivot assumptions needed to stay on a chosen fixed-rank symbolic stratum when a more
complicated Kossakowski form can change rank. Onset information and any method-specific branch data
belong to [`factorization`](@ref), not to the common physical interface.

The conceptual distinction matters in larger symbolic models: positivity asks whether a retained
Hermitian direction is physically admissible, whereas regularity asks whether the algebraic rank
pattern assumed by the current factorization remains valid.

See [Positive completion](@ref) for the API guide and
[High-frequency expansion](../theory/high_frequency_expansion.md) for the active/dark and onset
theory behind these examples.
