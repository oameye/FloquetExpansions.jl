struct OpenLegHoriDepritOrder3{P}
  generator1::P
  generator2::P
  generator3::P
  effective1::P
  effective2::P
  effective3::P
end

function openleg_hori_deprit_order3(A1::P, A2::P; product) where {P<:OpenLegPeriodic}
  lower = openleg_hori_deprit_order2(A1, A2; product)
  G1 = lower.generator1
  G2 = lower.generator2
  derivative_G1 = openleg_derivative(G1)
  derivative_G2 = openleg_derivative(G2)

  F3 =
    -openleg_commutator(G1, A2, product) - openleg_commutator(G2, A1, product) +
    (1 // 2) * openleg_commutator(G1, openleg_commutator(G1, A1, product), product) +
    (1 // 2) * openleg_commutator(G1, derivative_G2, product) +
    (1 // 2) * openleg_commutator(G2, derivative_G1, product) -
    (1 // 6) *
    openleg_commutator(G1, openleg_commutator(G1, derivative_G1, product), product)

  B3 = openleg_project(F3)
  G3 = openleg_q_inverse(F3)
  return OpenLegHoriDepritOrder3(G1, G2, G3, lower.effective1, lower.effective2, B3)
end
