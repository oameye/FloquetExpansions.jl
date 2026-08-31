using Documenter: DocMeta, doctest
using FloquetExpansions

# Documenter evaluates manual-page metadata in the worker's Main module.
Core.eval(Main, :(import FloquetExpansions))

DocMeta.setdocmeta!(
  FloquetExpansions, :DocTestSetup, :(using FloquetExpansions); recursive=true
)

doctest(FloquetExpansions)
