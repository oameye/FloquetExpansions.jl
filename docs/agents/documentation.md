# Documentation

Authority for how user-facing documentation is written: the pages under `docs/src/` and the docstrings of the public interface. `STANDARDS.md` routes here. The build, doctest, and generated-example workflow belongs to [`development.md`](development.md).

## Where each fact lives

Each fact lives in exactly one place, chosen by the question it answers. The other places link to it rather than restate it.

| Place | Answers | Holds |
| --- | --- | --- |
| Docstring | What does this name do? | What it selects, returns, or computes; its guarantees and caveats; one doctest |
| Manual page | Which name do I need, in what order, and how do I choose? | Short connective prose between `@docs` blocks |
| Theory page | Why does it work? | Conventions, defining equations, mechanism, literature |
| Example | What does real use look like? | A worked model, generated from `examples/*.jl` |

**Keep implementation out of user-facing documentation.** How the package organizes a computation, such as its internal representations, compiled plans, or performance strategy, belongs in [`architecture.md`](architecture.md), not in a docstring, manual page, or theory page.

## Docstrings

Every exported name and every name in the qualified expert API marked with `@public` has a docstring. A docstring is read alone at the REPL through `?name`, so it is **self-contained**: a reader who never opens the manual can use the name correctly. Self-contained means usable, not derived. Equations and mechanism belong to the theory page, reached through one link such as ``[High-frequency expansion](@ref high-frequency-expansion-theory)``. A guarantee or caveat about the name itself stays in its docstring, even when a manual reader would also want it.

## Manual pages

A page under `docs/src/manual/` walks through one stage of the user workflow, and it is **short**: usually one or two sentences before each `@docs` block, sometimes none. Its prose carries only what makes sense between docstrings: the order in which a user meets the names, comparisons between them, and how to choose. Section headings name user tasks, such as "Choosing a gauge" or "Reading the result". Two tests hold a page to this:

- **Paste test.** A sentence that would still make sense pasted into a docstring belongs in that docstring.
- **Deletion test.** With the `@docs` blocks removed, the remaining prose reads as a short list of decisions, not as a summary of the API.

## Theory pages

A page under `docs/src/theory/` is a **brief introduction** to the physics and mathematics behind a construction: conventions, defining equations, the mechanism, and the literature. It stops once the reader can see what the package computes and why, and leaves full derivations to the cited references. It uses the package conventions, such as the Fourier sign and the order indexing, and says where the literature differs. It holds no `@docs` blocks and describes no API; it links to the manual page for usage.

## Linking and citations

Give a new page an `@id` anchor so links to it survive a title change, and register it in `docs/pages.jl`. Cite with `[Key](@cite)` and add the entry to `docs/src/refs.bib`, taking its metadata from the published record: authors, journal, volume, pages, year, and DOI.
