CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using FloquetExpansions
using Documenter
using DocumenterCitations
using DocumenterCodeBlocks
using DocumenterInterLinks
using DocumenterLandingPage

using Plots
default(; fmt=:png)
# Gotta set this environment variable when using the GR run-time on CI machines.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

include("pages.jl")

bib = CitationBibliography("src/refs.bib")
links = InterLinks(
  "Julia" => "https://docs.julialang.org/en/v1/",
  "Documenter" => "https://documenter.juliadocs.org/stable/",
)

# The README.md file is used index (home) page of the documentation.
if CI
  include("make_md_examples.jl")
else
  nothing
end
# ^ when using LiveServer, this will generate a loop

DocMeta.setdocmeta!(
  FloquetExpansions,
  :DocTestSetup,
  :(using FloquetExpansions; using LinearAlgebra: ishermitian);
  recursive=true,
)

makedocs(;
  sitename="FloquetExpansions.jl",
  authors="Orjan Ameye",
  modules=[FloquetExpansions],
  format=Documenter.HTML(; canonical="https://oameye.github.io/FloquetExpansions.jl"),
  pages=pages,
  plugins=[bib, CodeBlocks(), LandingPage(), links],
  clean=true,
  linkcheck=false,
  draft=false,#,(!CI),
  doctest=false, # run in test suite
  checkdocs=:exports,
)

if CI
  deploydocs(;
    repo="github.com/oameye/FloquetExpansions.jl",
    devbranch="main",
    target="build",
    branch="gh-pages",
    push_preview=true,
  )
end
