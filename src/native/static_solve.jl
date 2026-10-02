struct NativeStaticSolution
  coordinates::Vector{Float64}
  newborn::Matrix{ComplexF64}
  correction::Matrix{ComplexF64}
  coefficient::Matrix{ComplexF64}
  dark_residual::Matrix{ComplexF64}
  canonical::Bool
  iterations::Int
end

function native_dark_section(
  residual::Matrix{ComplexF64},
  known::Matrix{ComplexF64},
  dark::Matrix{ComplexF64},
  gauge_images::Vector{Matrix{ComplexF64}},
  tol::Float64,
)
  delta = hermitian_part(adjoint(dark) * (residual - known) * dark)
  Phi = zeros(Float64, size(dark, 2)^2, length(gauge_images))
  for (column, image) in pairs(gauge_images)
    Phi[:, column] = gram_hermitian_coordinates(
      hermitian_part(adjoint(dark) * image * dark)
    )
  end
  section = positive_affine_section(delta, Phi; rtol=1e-10)
  form = Matrix{ComplexF64}(hermitian_part(section.form))
  newborn = dark * positive_gram_factor(form; rtol=tol)
  coordinates = Vector{Float64}(section.coordinates)
  return (coordinates, newborn, form, section.canonical, section.iterations)
end

function native_static_solve(
  residual::AbstractMatrix{<:Number},
  known::AbstractMatrix{<:Number},
  active::AbstractMatrix{<:Number},
  gauge_images::AbstractVector{<:AbstractMatrix{<:Number}},
  tol::Real,
)
  dimension = size(residual, 1)
  size(residual) == size(known) == (dimension, dimension) || throw(
    DimensionMismatch("residual and known Gram coefficients must be square and equal")
  )
  size(active, 1) == dimension ||
    throw(DimensionMismatch("active amplitudes must match the Kossakowski dimension"))
  V = Matrix{ComplexF64}(residual)
  K = Matrix{ComplexF64}(known)
  B = Matrix{ComplexF64}(active)
  images = Matrix{ComplexF64}[Matrix{ComplexF64}(image) for image in gauge_images]
  tolerance = Float64(tol)

  dark = Matrix{ComplexF64}(gram_active_frame(B; rtol=tolerance).dark)
  coordinates, newborn, dark_residual, canonical, iterations = if size(dark, 2) == 0
    (
      zeros(Float64, length(images)),
      zeros(ComplexF64, dimension, 0),
      zeros(ComplexF64, 0, 0),
      true,
      0,
    )
  else
    native_dark_section(V, K, dark, images, tolerance)
  end

  solved = copy(V)
  for (coordinate, image) in zip(coordinates, images)
    solved += coordinate * image
  end
  born = newborn * adjoint(newborn)
  correction = gram_tangent_lift(B, hermitian_part(solved - K - born); rtol=tolerance)
  coefficient = K + B * adjoint(correction) + correction * adjoint(B) + born
  LinearAlgebra.norm(solved - coefficient) <=
  tolerance * max(1.0, LinearAlgebra.norm(coefficient)) || throw(
    ArgumentError("native static slot failed to reconstruct its Kossakowski coefficient")
  )
  return NativeStaticSolution(
    coordinates, newborn, correction, coefficient, dark_residual, canonical, iterations
  )
end
