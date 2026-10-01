abstract type NativeRepresentation end

abstract type HomologicalInverse end

function native_kossakowski end
function native_hamiltonian end
function native_gksl end
function native_gauge_directions end

struct NoHomologicalInverse <: HomologicalInverse end

function regular_gauge(
  ::NoHomologicalInverse, ::NativeRepresentation, L0, residual, known, dark
)
  return zero(L0)
end

function singular_gauge_directions(
  ::NoHomologicalInverse, representation::NativeRepresentation
)
  return native_gauge_directions(representation)
end

struct NativeStaticStep{X,O}
  E::X
  S::X
  H::O
  solution::NativeStaticSolution
end

native_commutator(A, B) = A * B - B * A

function native_static_step(
  representation::NativeRepresentation,
  inverse::HomologicalInverse,
  L0,
  residual,
  known::AbstractMatrix{<:Number},
  active::AbstractMatrix{<:Number},
  tol::Real,
)
  dark = gram_active_frame(Matrix{ComplexF64}(active); rtol=tol).dark
  regular = regular_gauge(inverse, representation, L0, residual, known, dark)
  singular_residual = residual + native_commutator(L0, regular)
  directions = singular_gauge_directions(inverse, representation)
  images = Matrix{ComplexF64}[
    native_kossakowski(representation, native_commutator(L0, G)) for G in directions
  ]
  solution = native_static_solve(
    native_kossakowski(representation, singular_residual), known, active, images, tol
  )

  S = regular
  for (coordinate, G) in zip(solution.coordinates, directions)
    S += coordinate * G
  end
  E = residual + native_commutator(L0, S)

  coefficient = solution.coefficient
  scale = max(1.0, LinearAlgebra.norm(coefficient))
  LinearAlgebra.norm(native_kossakowski(representation, E) - coefficient) <= tol * scale ||
    throw(
      ArgumentError("native static step failed to reconstruct its Kossakowski coefficient")
    )
  H = native_hamiltonian(representation, E)
  reconstruction = native_gksl(representation, H, coefficient)
  LinearAlgebra.norm(E - reconstruction) <= 20 * tol * max(1.0, LinearAlgebra.norm(E)) ||
    throw(
      ArgumentError("native static step failed its Hamiltonian/Kossakowski reconstruction")
    )
  return NativeStaticStep(E, S, H, solution)
end
