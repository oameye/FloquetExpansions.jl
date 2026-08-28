JULIA:=julia

# `test`, `docs` and `examples` are also directory names; without this Make reports
# "'test' is up to date" and runs nothing.
.PHONY: default setup format servedocs test jet docs bench all help

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
	@echo " - make all: run every commands in the above order"

