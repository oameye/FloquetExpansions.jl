JULIA:=julia

# `test`, `docs` and `examples` are also directory names; without this Make reports
# "'test' is up to date" and runs nothing.
.PHONY: default setup format servedocs test jet docs bench ratchet ratchet-scorecard ratchet-candidates ratchet-coverage ratchet-boxes ratchet-undocumented ratchet-lsp-report ratchet-refresh all help

default: help

setup:
	${JULIA} -e 'import Pkg; Pkg.add(["Changelog", "LiveServer"])'
	${JULIA} -e 'using Pkg; Pkg.Apps.add("JuliaFormatter")'


format: ## Format all Julia files with JuliaFormatter
	jlfmt --threads=6 -- --inplace -v ./

servedocs:
	${JULIA} --project=docs -e 'using LiveServer; LiveServer.servedocs(skip_files=[joinpath("docs", "src", "changelog.md")])'

test:
	${JULIA} --project -e 'using Pkg; Pkg.resolve(); Pkg.test()'

jet:
	${JULIA} --project=test -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=test test/quality/JET.jl

docs:
	${JULIA} --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=docs docs/make.jl

bench:
	${JULIA} --project=benchmark -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=benchmark benchmark/runbenchmarks.jl

# The code-quality ratchet. Each number only has to stop getting worse, so these
# are green on a clean tree and go red when a file regresses.
#
# Which gates run is [metrics].run in code_ratchet/rulings.toml, not a list
# here: two places to write the set down is one place too many. `all` runs them
# cheapest first, so the second-long gates fail before the minutes-long one.
ratchet:
	${JULIA} --project=code_ratchet -e 'using Pkg; Pkg.instantiate()'
	${JULIA} --project=code_ratchet -e 'using JET, CodeRatchet; exit(CodeRatchet.main())' all check

# Where the recorded debt is, rather than what regressed. Reads the baselines,
# so it is instant and gates nothing.
ratchet-scorecard:
	${JULIA} --project=code_ratchet -e 'using JET, CodeRatchet; exit(CodeRatchet.main())' all scorecard

# What to fix, named rather than counted. None of these gate.
ratchet-candidates:
	${JULIA} --project=code_ratchet -e 'using CodeRatchet; exit(CodeRatchet.main())' complexity candidates

ratchet-boxes:
	${JULIA} --project=code_ratchet -e 'using CodeRatchet; exit(CodeRatchet.main())' boxes methods

ratchet-undocumented:
	${JULIA} --project=code_ratchet -e 'using CodeRatchet; exit(CodeRatchet.main())' docs undocumented

ratchet-lsp-report:
	${JULIA} --project=code_ratchet -e 'using CodeRatchet; exit(CodeRatchet.main())' lsp report

# Coverage needs an lcov tracefile, which the test job produces. Point
# COVERAGE_LCOV at one, or drop it beside this Makefile as lcov.info. It is not
# in [metrics].run for that reason.
ratchet-coverage:
	${JULIA} --project=code_ratchet -e 'using CodeRatchet; exit(CodeRatchet.main())' coverage check

# Rewrites every baseline. Refuses a change unless --accept-change is passed
# through, so recording a worse number stays a deliberate act.
ratchet-refresh:
	${JULIA} --project=code_ratchet -e 'using JET, CodeRatchet; exit(CodeRatchet.main())' all refresh

all: setup format test docs

help:
	@echo "The following make commands are available:"
	@echo " - make setup: install the dependencies for make command"
	@echo " - make format: format codes with JuliaFormatter"
	@echo " - make test: run the tests"
	@echo " - make jet: run JET static analysis"
	@echo " - make docs: instantiate and build the documentation"
	@echo " - make servedocs: serve the documentation locally"
	@echo " - make bench: run the benchmarks"
	@echo " - make ratchet: run every configured gate, cheapest first"
	@echo " - make ratchet-scorecard: show where the recorded debt is"
	@echo " - make ratchet-candidates: rank definitions above threshold"
	@echo " - make ratchet-boxes: name every boxed closure capture"
	@echo " - make ratchet-undocumented: name every public name owing a docstring"
	@echo " - make ratchet-lsp-report: list every undismissed JETLS diagnostic"
	@echo " - make ratchet-coverage: run the coverage ratchet (needs lcov.info)"
	@echo " - make ratchet-refresh: rewrite every ratchet baseline"
	@echo " - make all: run every commands in the above order"

