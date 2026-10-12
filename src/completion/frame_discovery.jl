frame_seeds(::NoProvenance) = SQA.QAdd[]
frame_seeds(provenance::MicroscopicProvenance) = provenance.frame_seeds

function automatic_dissipative_frame(expansion::FloquetExpansion)
  candidates = copy(frame_seeds(getfield(expansion, :provenance)))
  for component in getfield(expansion, :effective_components)
    append!(candidates, SQA.directions(SQA.support_frame(component)))
  end
  return SQA.independent_frame(candidates)
end
