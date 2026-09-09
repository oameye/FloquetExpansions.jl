# Julia Style

Authority for how source in `src/` is written: signatures, fields, imports, formatting, comments, and names. `STANDARDS.md` routes here and names the gate, if any, for each rule below.

Scope: `src/`. Test code follows `development.md`.

## Function signatures

- **Constrain a parameter when the constraint expresses a semantic invariant or dispatch boundary.** Do not make a signature artificially concrete merely to silence inference tooling. The expansion engine is intentionally generic over the `PeriodicGenerator` component algebra, while completion kernels use stronger constraints where the representation requires them.
- **Write keyword arguments after an explicit `;`**, at the call site as well as the definition:

  ```julia
  FockSpace(; name = :a)   # yes
  FockSpace(name = :a)     # no
  ```

- **Use keyword shorthand when forwarding same-named locals.** Write `f(args...; kwarg1, kwarg2)` rather than `f(args...; kwarg1 = kwarg1, kwarg2 = kwarg2)`.
- **Leave a higher-order function's parameter untyped by default.** Use `f::F where {F}` when specialization is deliberately required, such as a small callback on a compiler-sensitive inner loop.

## Type system

**Keep struct storage concrete.** Fields may be concrete types or type parameters that become concrete for each instance; do not store values behind `Any` or an abstract field type when the representation can encode the concrete type. Gated by `test/quality/CheckConcreteStructs.jl`.

## Imports

**Import explicitly:** `using X: func1, func2`, or `import X`. Gated by `test/quality/ExplicitImports.jl`, which also rejects stale explicit imports, owner-incorrect imports or qualified accesses, and self-qualified accesses. `Base`, `Core`, and `SecondQuantizedAlgebra` are skipped only for the implicit-import check configured by that test.

## Formatting

JuliaFormatter owns formatting, configured by `.JuliaFormatter.toml`: `blue` style, indent 2. Run `make format`. `.github/workflows/Format.yml` runs `jlfmt --check --verbose .` over the repository.

## Comments and docstrings

- **Comment the non-obvious why, not the visible what.** Keep comments compact and local to the invariant they explain.
- **Docstring the public interface.** This includes exported names and the intentionally qualified expert API marked with `@public`. Ordinary internal helpers should normally be explained by their names, types, and the architecture/ADR that owns any non-obvious invariant rather than by duplicating design prose in local docstrings.

## Naming

**Name private internals like ordinary Julia names.** A leading underscore (`_helper`, `_util.jl`) is not the repository's convention for privacy; module visibility and the explicit export/`@public` surface define the API.

Take a domain concept's name from `CONTEXT.md`. Where the concept is not there yet, add it in the same change.
