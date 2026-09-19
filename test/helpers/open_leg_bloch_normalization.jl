function openleg_bloch_static_factor2(bloch::OpenLegBlochReference; product)
  X1 = bloch.wave[1]
  return (1 // 2) * openleg_project(openleg_product(X1, X1, product))
end

function openleg_bloch_canonical_effective3(bloch::OpenLegBlochReference; product)
  length(bloch.effective) >= 3 ||
    throw(ArgumentError("third-order Bloch data are required"))
  static2 = openleg_bloch_static_factor2(bloch; product)
  correction = openleg_commutator(bloch.effective[1], static2, product)
  return bloch.effective[3] + correction
end
