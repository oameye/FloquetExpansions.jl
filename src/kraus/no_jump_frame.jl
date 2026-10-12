struct NoJumpFrame{X}
  effective::Vector{X}
  kick::Vector{Harmonics{X}}
  inverse::Vector{Harmonics{X}}
end

function kick_inverse(K::Vector{Harmonics{X}}, identity::X, N::Int) where {X}
  A = kick_add(kick_identity(identity, N), K, N, 1, -1)
  result = kick_identity(identity, N)
  power = kick_identity(identity, N)
  for _ in 1:N
    power = kick_mul(power, A, N)
    result = kick_add(result, power, N, 1, 1)
  end
  return result
end

function frame_kicks(::HoriDeprit, G::Vector{<:Harmonics}, identity, N::Int)
  return kick_exp(G, identity, N), kick_exp(kick_scale(G, N, -1), identity, N)
end

function frame_kicks(::BlochFeshbach, Y::Vector{<:Harmonics}, identity, N::Int)
  return Y, kick_inverse(Y, identity, N)
end

function no_jump_frame(
  algorithm::Union{BlochFeshbach,HoriDeprit}, Q::Harmonics{X}, N::Int
) where {X}
  effective, _, kick = floquet_recurrence(algorithm, Q, N, VanVleckSlot())
  return NoJumpFrame{X}(effective, frame_kicks(algorithm, kick, one(Q[0]), N)...)
end

function endpoint_value(series::Vector{Harmonics{X}}, zero_element::X) where {X}
  return X[sum(values(harmonics); init=zero_element) for harmonics in series]
end

function transported_jump(frame::NoJumpFrame{X}, J::Harmonics{X}, N::Int) where {X}
  series = kick_zero(X, N)
  series[1] = J
  return kick_mul(kick_mul(frame.inverse, series, N), frame.kick, N)
end

function series_product(A::Vector{X}, B::Vector{X}, N::Int, zero_element::X) where {X}
  result = X[zero_element for _ in 0:N]
  for (i, a) in pairs(A), (j, b) in pairs(B)
    i + j - 2 <= N || continue
    result[i + j - 1] += a * b
  end
  return result
end

function effective_powers(frame::NoJumpFrame{X}, N::Int) where {X}
  identity = one(frame.effective[1])
  zero_element = zero(identity)
  power = X[n == 0 ? identity : zero_element for n in 0:N]
  powers = [power]
  for p in 1:N
    power = series_product(power, frame.effective, N - p, zero_element)
    push!(powers, power)
  end
  return powers
end
