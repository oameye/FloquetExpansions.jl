# Domain Docs

How the engineering skills should consume this repo's domain documentation.

## Before exploring, read these

- **`CONTEXT.md`** at the repo root.
- **`docs/adr/`**: read ADRs that touch the area being explored.

If these files don't exist, proceed silently. The `/domain-modeling` skill creates them lazily when terms or decisions are resolved.

## File structure

Single-context repository:

/
├── CONTEXT.md
├── docs/adr/
└── src/

## Use the glossary's vocabulary

When naming a domain concept in an issue, proposal, or test, use the term defined in `CONTEXT.md`. If the concept is not yet defined, note it for `/domain-modeling`.

## Every ADR names its own gate

An ADR ends with a `## Gate` section listing the testsets that hold its decision, by file and
testset name, or stating that no gate holds it. The gate belongs with the decision: a decision
whose text lives in one file and whose check is named in another is one rule in two places, and
the two drift. `STANDARDS.md` therefore routes no ADR.

Adding an ADR means adding that section. Changing which tests hold a decision means updating it.

## Flag ADR conflicts

If output contradicts an existing ADR, surface it explicitly rather than silently overriding it.
