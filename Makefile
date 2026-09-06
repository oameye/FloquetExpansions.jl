JULIA:=julia

# `test`, `docs` and `examples` are also directory names; without this Make reports
# "'test' is up to date" and runs nothing.
.PHONY: default setup format servedocs test jet docs bench ratchet ratchet-coverage ratchet-candidates ratchet-refresh all help

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
# are green on a clean tree and go red when a file regresses. `candidates` ranks
# work above the rulings.toml thresholds and never gates.
ratchet:
	${JULIA} --project=code_ratchet -e 'using Pkg; Pkg.instantiate()'
	${JULIA} --project=code_ratchet -e 'using CodeRatchet; exit(CodeRatchet.main())' complexity check
	${JULIA} --project=code_ratchet -e 'using JET, CodeRatchet; exit(CodeRatchet.main())' jet check

# Coverage needs an lcov tracefile, which the test job produces. Point
# COVERAGE_LCOV at one, or drop it beside this Makefile as lcov.info.
ratchet-coverage:
	${JULIA} --project=code_ratchet -e 'using CodeRatchet; exit(CodeRatchet.main())' coverage check

ratchet-candidates:
	${JULIA} --project=code_ratchet -e 'using CodeRatchet; exit(CodeRatchet.main())' complexity candidates

# Rewrites the baselines. Refuses a rise unless you pass --accept-rise through,
# so recording a worse number stays a deliberate act.
ratchet-refresh:
	${JULIA} --project=code_ratchet -e 'using CodeRatchet; exit(CodeRatchet.main())' complexity refresh
	${JULIA} --project=code_ratchet -e 'using JET, CodeRatchet; exit(CodeRatchet.main())' jet refresh

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
	@echo " - make ratchet: run the code-quality ratchet (complexity + JET)"
	@echo " - make ratchet-coverage: run the coverage ratchet (needs lcov.info)"
	@echo " - make ratchet-candidates: rank definitions above threshold"
	@echo " - make ratchet-refresh: rewrite the ratchet baselines"
	@echo " - make all: run every commands in the above order"

