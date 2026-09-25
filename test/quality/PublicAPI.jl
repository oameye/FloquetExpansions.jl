using Test
using FloquetExpansions

const PUBLIC_API_TEST_ROOT = normpath(joinpath(@__DIR__, ".."))
const PUBLIC_API_MODULE_FILE = normpath(
  joinpath(PUBLIC_API_TEST_ROOT, "..", "src", "FloquetExpansions.jl")
)
const PUBLIC_API_SELF = normpath(@__FILE__)

function declared_public_api()
  public_names = Set{Symbol}(names(FloquetExpansions; all=false, imported=false))
  source = read(PUBLIC_API_MODULE_FILE, String)
  for block in eachmatch(r"(?ms)@public\s+(.*?)(?=\n\n|\z)", source)
    for name in eachmatch(r"[A-Za-z_][A-Za-z0-9_!]*", only(block.captures))
      push!(public_names, Symbol(name.match))
    end
  end
  return public_names
end

function test_julia_files()
  files = String[]
  for (root, _, names) in walkdir(PUBLIC_API_TEST_ROOT)
    for name in names
      endswith(name, ".jl") || continue
      push!(files, normpath(joinpath(root, name)))
    end
  end
  return files
end

function imported_floquet_names(source::String)
  imported = Symbol[]
  pattern = r"(?ms)\b(?:using|import)\s+FloquetExpansions\s*:\s*(.*?)(?=\n\S|\z)"
  for block in eachmatch(pattern, source)
    for name in eachmatch(r"[A-Za-z_][A-Za-z0-9_!]*", only(block.captures))
      push!(imported, Symbol(name.match))
    end
  end
  return imported
end

@testset "tests use only the FloquetExpansions public API" begin
  allowed = declared_public_api()
  violations = String[]
  qualified_pattern = r"\bFloquetExpansions\.([A-Za-z_][A-Za-z0-9_!]*)"
  alias_pattern = r"(?m)^\s*(?:const\s+)?[A-Za-z_][A-Za-z0-9_]*\s*=\s*FloquetExpansions\s*$"

  for path in test_julia_files()
    source = read(path, String)
    relative = relpath(path, PUBLIC_API_TEST_ROOT)

    for name in imported_floquet_names(source)
      if !(name in allowed)
        push!(violations, "$relative imports private FloquetExpansions.$name")
      end
    end

    for qualified in eachmatch(qualified_pattern, source)
      name = Symbol(only(qualified.captures))
      if !(name in allowed)
        push!(violations, "$relative reaches private FloquetExpansions.$name")
      end
    end

    if occursin(alias_pattern, source)
      push!(
        violations, "$relative aliases the FloquetExpansions module and bypasses API checks"
      )
    end

    path == PUBLIC_API_SELF && continue
    if occursin("getfield(", source)
      push!(violations, "$relative uses getfield; tests must use public accessors")
    end
  end

  isempty(violations) || @error "private FloquetExpansions test access" violations
  @test isempty(violations)
end
