using Test
using LinearAlgebra: I, det, kron

const GramOracleExact = Complex{Rational{Int}}
const gram_oracle_im = GramOracleExact(0 // 1, 1 // 1)

function gram_oracle_identity(dimension::Int)
  return Matrix{GramOracleExact}(I, dimension, dimension)
end

function gram_oracle_metric(gram, amplitudes)
  result = zero(first(amplitudes))
  for alpha in eachindex(amplitudes), beta in eachindex(amplitudes)
    result += gram[alpha, beta] * adjoint(amplitudes[alpha]) * amplitudes[beta]
  end
  return result
end

function gram_oracle_channel(gram, amplitudes)
  dimension = size(first(amplitudes), 1)
  result = zeros(GramOracleExact, dimension^2, dimension^2)
  for alpha in eachindex(amplitudes), beta in eachindex(amplitudes)
    result += gram[beta, alpha] * kron(conj.(amplitudes[beta]), amplitudes[alpha])
  end
  return result
end

function gram_oracle_kraus(factor, amplitudes)
  rows = Matrix{GramOracleExact}[]
  for row in axes(factor, 1)
    kraus = zero(first(amplitudes))
    for alpha in eachindex(amplitudes)
      kraus += factor[row, alpha] * amplitudes[alpha]
    end
    push!(rows, kraus)
  end
  return rows
end

function gram_oracle_kraus_metric(rows)
  result = zero(first(rows))
  for row in rows
    result += adjoint(row) * row
  end
  return result
end

function gram_oracle_kraus_channel(rows)
  dimension = size(first(rows), 1)
  result = zeros(GramOracleExact, dimension^2, dimension^2)
  for row in rows
    result += kron(conj.(row), row)
  end
  return result
end

@testset "nonorthogonal Gram channel equals an independent explicit Kraus factorization" begin
  factor = GramOracleExact[1 gram_oracle_im; 1 0]
  gram = adjoint(factor) * factor
  @test gram == GramOracleExact[2 gram_oracle_im; -gram_oracle_im 1]

  amplitudes = [GramOracleExact[1 1; 0 -1], GramOracleExact[0 1 + gram_oracle_im; 1 -1]]
  kraus = gram_oracle_kraus(factor, amplitudes)

  @test gram_oracle_metric(gram, amplitudes) == gram_oracle_kraus_metric(kraus)
  @test gram_oracle_channel(gram, amplitudes) == gram_oracle_kraus_channel(kraus)

  # `mixed_qr_gram_polynomial(left,right)` uses the reverse-bra convention
  # G_RL[alpha,beta] = <e_beta|e_alpha> = G[beta,alpha].
  gram_right_left = transpose(gram)
  metric_from_right_left = zero(first(amplitudes))
  channel_from_right_left = zeros(GramOracleExact, 4, 4)
  for alpha in eachindex(amplitudes), beta in eachindex(amplitudes)
    metric_from_right_left +=
      conj(gram_right_left[alpha, beta]) * adjoint(amplitudes[alpha]) * amplitudes[beta]
    channel_from_right_left +=
      gram_right_left[alpha, beta] * kron(conj.(amplitudes[beta]), amplitudes[alpha])
  end
  @test metric_from_right_left == gram_oracle_kraus_metric(kraus)
  @test channel_from_right_left == gram_oracle_kraus_channel(kraus)
end

@testset "singular physical Gram basis normalizes exactly without inverting its Gram matrix" begin
  factor = reshape(GramOracleExact[1, gram_oracle_im], 1, 2)
  gram = adjoint(factor) * factor
  @test gram == GramOracleExact[1 gram_oracle_im; -gram_oracle_im 1]
  @test iszero(det(gram))

  identity_component = gram_oracle_identity(2)
  amplitudes = [identity_component, -gram_oracle_im * identity_component]
  metric = gram_oracle_metric(gram, amplitudes)
  @test metric == 4 * identity_component

  normalization = (1 // 2) * identity_component
  normalized_amplitudes = [amplitude * normalization for amplitude in amplitudes]
  @test gram_oracle_metric(gram, normalized_amplitudes) == identity_component

  kraus = gram_oracle_kraus(factor, amplitudes)
  normalized_kraus = gram_oracle_kraus(factor, normalized_amplitudes)
  @test kraus == [2 * identity_component]
  @test normalized_kraus == [identity_component]
  @test gram_oracle_kraus_metric(normalized_kraus) == identity_component
  @test gram_oracle_kraus_channel(normalized_kraus) == Matrix{GramOracleExact}(I, 4, 4)
end
