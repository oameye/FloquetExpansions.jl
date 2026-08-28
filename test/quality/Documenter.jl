using Documenter: DocMeta, doctest
using FloquetExpansions

DocMeta.setdocmeta!(
  FloquetExpansions,
  :DocTestSetup,
  :(using FloquetExpansions),
  recursive=true,
)

doctest(FloquetExpansions)
