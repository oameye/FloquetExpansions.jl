# Julia Style

Authority for how source in `src/` is written: signatures, types, imports, comments, names. `STANDARDS.md` routes here and names the gate for each rule below.

Scope: `src/**/*.jl`. Test code follows `development.md`.

## Function signatures

- **Use the most restrictive signature type possible.** Tight type declarations let JET catch errors. Starting loose while prototyping is fine; committed code carries specific types.
- **Write keyword arguments after an explicit `;`**, at the call site as well as the definition:

  ```julia
  FockSpace(; name = :a)   # yes
  FockSpace(name = :a)     # no
  ```

- **Use keyword shorthand when forwarding same-named locals.** Write `f(args...; kwarg1, kwarg2)` rather than `f(args...; kwarg1 = kwarg1, kwarg2 = kwarg2)`.
- **Leave a higher-order function's parameter untyped.** Julia's pass-through heuristic skips specialization for a `Function`-annotated argument that is not called directly, and the same heuristic applies to an untyped `f`. Reach for `f::F where {F}` only to force specialization, such as a nested closure. The difference is usually negligible, so keep the signature simple.

## Type system

**Give every struct field a concrete type.** Gated by `test/quality/CheckConcreteStructs.jl`.

## Imports

**Import explicitly:** `using X: func1, func2`, or `import X`. Gated by `test/quality/ExplicitImports.jl`, which also rejects a stale explicit import, an import that bypasses the owning module, and a self-qualified access. `Base`, `Core` and `SecondQuantizedAlgebra` are skipped for the implicit-import check.

## Formatting

JuliaFormatter owns formatting, configured by `.JuliaFormatter.toml`: `blue` style, indent 2. Run `make format`. Gated by `.github/workflows/Format.yml`, which runs `jlfmt --check --verbose .` over the whole repository.

## Comments

- **Comment the non-obvious *why*, in a line or two.** The default is no comment, because the code states what it does. A comment restating the code is the one to delete.
- **Docstring the public interface.** An internal helper's name and signature carry it; a docstring there is load with no reader.

## Naming

**Name private internals like anything else.** Julia's module system handles visibility, so a leading underscore (`_helper`, `_util.jl`) buys nothing here.

Take a domain concept's name from `CONTEXT.md`. Where the concept is not there yet, add it in the same change.
