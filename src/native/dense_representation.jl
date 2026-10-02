struct DenseLiouvilleRepresentation <: NativeRepresentation
  d::Int
  basis::Vector{Matrix{ComplexF64}}
  gauge::Vector{Matrix{ComplexF64}}
end

function DenseLiouvilleRepresentation(d::Int)
  d >= 2 || throw(ArgumentError("a dense Liouville representation needs dimension d >= 2"))
  basis = dense_traceless_basis(d)
  return DenseLiouvilleRepresentation(d, basis, dense_gauge_algebra(d, basis))
end

dense_identity(d::Int) = Matrix{ComplexF64}(LinearAlgebra.I, d, d)
dense_left(A::AbstractMatrix) = kron(dense_identity(size(A, 1)), A)
dense_right(B::AbstractMatrix) = kron(transpose(B), dense_identity(size(B, 1)))
dense_sandwich(X::AbstractMatrix, Y::AbstractMatrix) = kron(conj(Y), X)

function dense_traceless_basis(d::Int)
  basis = Matrix{ComplexF64}[]
  for j in 1:d, k in (j + 1):d
    X = zeros(ComplexF64, d, d)
    X[j, k] = 1
    X[k, j] = 1
    push!(basis, X / sqrt(2))

    Y = zeros(ComplexF64, d, d)
    Y[j, k] = -im
    Y[k, j] = im
    push!(basis, Y / sqrt(2))
  end
  for l in 1:(d - 1)
    Z = zeros(ComplexF64, d, d)
    for j in 1:l
      Z[j, j] = 1
    end
    Z[l + 1, l + 1] = -l
    push!(basis, Z / sqrt(l * (l + 1)))
  end
  return basis
end

function dense_full_basis(d::Int, basis::Vector{Matrix{ComplexF64}})
  return vcat([dense_identity(d) / sqrt(d)], basis)
end

function dense_chi_matrix(d::Int, basis::Vector{Matrix{ComplexF64}}, L::AbstractMatrix)
  full = dense_full_basis(d, basis)
  n = d^2
  chi = zeros(ComplexF64, n, n)
  for i in 1:n, j in 1:n
    chi[i, j] = LinearAlgebra.dot(vec(dense_sandwich(full[i], full[j])), vec(L))
  end
  return chi
end

function dense_superoperator(full::Vector{Matrix{ComplexF64}}, chi::AbstractMatrix)
  n = length(full)
  result = zeros(ComplexF64, n, n)
  for i in 1:n, j in 1:n
    result += chi[i, j] * dense_sandwich(full[i], full[j])
  end
  return result
end

function dense_trace_operator(full::Vector{Matrix{ComplexF64}}, chi::AbstractMatrix)
  d = size(full[1], 1)
  result = zeros(ComplexF64, d, d)
  for i in eachindex(full), j in eachindex(full)
    result += chi[i, j] * adjoint(full[j]) * full[i]
  end
  return hermitian_part(result)
end

function dense_gauge_algebra(d::Int, basis::Vector{Matrix{ComplexF64}})
  full = dense_full_basis(d, basis)
  n = d^2
  columns = Vector{Vector{Float64}}()
  for k in 1:(n ^ 2)
    unit = zeros(Float64, n^2)
    unit[k] = 1
    chi = gram_hermitian_from_coordinates(unit, n)
    push!(columns, gram_hermitian_coordinates(dense_trace_operator(full, chi)))
  end
  kernel = LinearAlgebra.nullspace(reduce(hcat, columns))
  return [
    dense_superoperator(full, gram_hermitian_from_coordinates(kernel[:, j], n)) for
    j in axes(kernel, 2)
  ]
end

function native_kossakowski(representation::DenseLiouvilleRepresentation, L::AbstractMatrix)
  return dense_chi_matrix(representation.d, representation.basis, L)[2:end, 2:end]
end

function native_hamiltonian(representation::DenseLiouvilleRepresentation, L::AbstractMatrix)
  d = representation.d
  chi = dense_chi_matrix(d, representation.basis, L)
  H = zeros(ComplexF64, d, d)
  for (i, F) in pairs(representation.basis)
    H += chi[i + 1, 1] / sqrt(d) * F
  end
  return hermitian_part((im / 2) * (H - adjoint(H)))
end

function native_gksl(
  representation::DenseLiouvilleRepresentation, H::AbstractMatrix, C::AbstractMatrix
)
  basis = representation.basis
  result = -im * (dense_left(H) - dense_right(H))
  for i in eachindex(basis), j in eachindex(basis)
    abs(C[i, j]) <= 1e-14 && continue
    F = adjoint(basis[j]) * basis[i]
    result +=
      C[i, j] *
      (dense_sandwich(basis[i], basis[j]) - dense_left(F) / 2 - dense_right(F) / 2)
  end
  return result
end

native_gauge_directions(representation::DenseLiouvilleRepresentation) = representation.gauge
