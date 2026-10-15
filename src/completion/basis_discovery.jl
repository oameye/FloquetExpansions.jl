basis_seeds(::NoProvenance) = SQA.QAdd[]
basis_seeds(provenance::MicroscopicProvenance) = provenance.basis_seeds

function automatic_jump_basis(expansion::FloquetExpansion)
  candidates = copy(basis_seeds(getfield(expansion, :provenance)))
  for component in getfield(expansion, :effective_components)
    append!(candidates, SQA.support_basis(component))
  end
  return SQA.independent_basis(candidates)
end
