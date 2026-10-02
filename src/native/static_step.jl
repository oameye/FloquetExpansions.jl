abstract type NativeRepresentation end

abstract type HomologicalInverse end

function native_kossakowski end
function native_hamiltonian end
function native_gksl end
function native_gauge_directions end

struct NoHomologicalInverse <: HomologicalInverse end

function regular_gauge(
  ::NoHomologicalInverse, ::NativeRepresentation, L0, residual, known, active, tol
)
  return zero(L0)
end

function singular_gauge_directions(
  ::NoHomologicalInverse, representation::NativeRepresentation
)
  return native_gauge_directions(representation)
end

struct NativeStaticStep{X,O,S,C,B}
  E::X
  S::X
  H::O
  solution::S
  correction::C
  births::B
end

native_commutator(A, B) = A * B - B * A

function native_matches(::NativeRepresentation, A, B, tol::Real)
  return LinearAlgebra.norm(A - B) <= tol * max(1.0, LinearAlgebra.norm(B))
end

function native_dark_target(::NativeRepresentation, residual, known, active, tol::Real)
  dark = Matrix{ComplexF64}(gram_active_frame(Matrix{ComplexF64}(active); rtol=tol).dark)
  size(dark, 2) == 0 && return zeros(ComplexF64, size(residual))
  delta = hermitian_part(adjoint(dark) * (residual - known) * dark)
  return dark * delta * adjoint(dark)
end

function native_slot_solve(
  ::NativeRepresentation, residual, known, active, weights, images, directions, tol::Real
)
  scale = LinearAlgebra.Diagonal(sqrt.(real.(weights)))
  solution = native_static_solve(residual, known, active * scale, images, tol)
  correction = solution.correction / scale
  births = ones(ComplexF64, size(solution.newborn, 2))
  return solution, correction, births
end

function native_static_step(
  representation::NativeRepresentation,
  inverse::HomologicalInverse,
  L0,
  residual,
  known::AbstractMatrix{<:Number},
  active::AbstractMatrix{<:Number},
  weights::AbstractVector{<:Number},
  tol::Real,
)
  C = native_kossakowski(representation, residual)
  regular = regular_gauge(inverse, representation, L0, C, known, active, tol)
  singular_residual = residual + native_commutator(L0, regular)
  directions = singular_gauge_directions(inverse, representation)
  images = [
    native_kossakowski(representation, native_commutator(L0, G)) for G in directions
  ]
  solution, correction, births = native_slot_solve(
    representation,
    native_kossakowski(representation, singular_residual),
    known,
    active,
    weights,
    images,
    directions,
    tol,
  )

  S = regular
  for (coordinate, G) in zip(solution.coordinates, directions)
    S += coordinate * G
  end
  E = residual + native_commutator(L0, S)

  coefficient = solution.coefficient
  native_matches(representation, native_kossakowski(representation, E), coefficient, tol) ||
    throw(
      ArgumentError("native static step failed to reconstruct its Kossakowski coefficient")
    )
  H = native_hamiltonian(representation, E)
  native_matches(
    representation, E, native_gksl(representation, H, coefficient), 20 * tol
  ) || throw(
    ArgumentError("native static step failed its Hamiltonian/Kossakowski reconstruction")
  )
  return NativeStaticStep(E, S, H, solution, correction, births)
end
