abstract type NativeRepresentation end

function native_kossakowski end
function native_hamiltonian end
function native_gksl end

struct NativeStaticStep{X,O,S}
  E::X
  S::X
  H::O
  solution::S
  virtual::Bool
end

native_commutator(A, B) = A * B - B * A

struct NativeGaugeFamily{D,I,M}
  directions::D
  images::I
  metric::M
end

native_gauge_metric(::NativeRepresentation, _) = zeros(Float64, 0, 0)

function native_gauge_family(representation::NativeRepresentation, directions, L0)
  images = [
    native_kossakowski(representation, native_commutator(L0, G)) for G in directions
  ]
  return NativeGaugeFamily(
    directions, images, native_gauge_metric(representation, directions)
  )
end

function native_static_step(
  representation::NativeRepresentation,
  family::NativeGaugeFamily,
  L0,
  residual,
  known::AbstractMatrix{<:Number},
  active::AbstractMatrix{<:Number},
  weights::AbstractVector{<:Number},
)
  directions = family.directions
  solution = native_slot_solve(
    representation,
    native_kossakowski(representation, residual),
    known,
    active,
    weights,
    family,
  )

  S = zero(L0)
  for (coordinate, G) in zip(solution.coordinates, directions)
    S += coordinate * G
  end
  E = residual + native_commutator(L0, S)

  coefficient = solution.coefficient
  native_matches(representation, native_kossakowski(representation, E), coefficient) ||
    throw(
      ArgumentError("native static step failed to reconstruct its Kossakowski coefficient")
    )
  H = native_hamiltonian(representation, E)
  native_matches(representation, E, native_gksl(representation, H, coefficient)) || throw(
    ArgumentError("native static step failed its Hamiltonian/Kossakowski reconstruction")
  )
  return NativeStaticStep(E, S, H, solution, false)
end

function native_static_step(
  representation::NativeRepresentation,
  families::AbstractVector{<:NativeGaugeFamily},
  L0,
  residual,
  known::AbstractMatrix{<:Number},
  active::AbstractMatrix{<:Number},
  weights::AbstractVector{<:Number},
)
  for family in families[1:(end - 1)]
    try
      return native_static_step(
        representation, family, L0, residual, known, active, weights
      )
    catch error
      error isa NativePositivityError || rethrow()
    end
  end
  return native_static_step(
    representation, families[end], L0, residual, known, active, weights
  )
end
