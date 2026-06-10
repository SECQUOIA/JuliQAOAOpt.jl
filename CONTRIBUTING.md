# Contributing

Thank you for contributing to JuliQAOAOpt.jl. This repository is a QUBODrivers
interface package for JuliQAOA, so changes should stay focused on the interface,
metadata, documentation, tests, and compatibility glue needed by this package.

For the generic mechanics of implementing a QUBODrivers sampler, use the
[QUBODrivers sampler setup guide](https://juliaqubo.github.io/QUBODrivers.jl/stable/manual/4-setup/).

## Supported Julia Versions

JuliQAOAOpt currently targets Julia 1.10 LTS and Julia 1.11 for angle finding.
The `JuliQAOA.find_angles_bh` path relies on upstream JuliQAOA and Enzyme support
that remains limited on Julia 1.12+. Keep compatibility and CI changes aligned
with `Project.toml` unless the upstream angle-search path changes.

## Development Setup

`JuliQAOA.jl` is not currently registered in Julia's General registry, so develop
it by URL before instantiating and testing this package:

```shell
julia --project=. -e 'import Pkg; Pkg.develop(Pkg.PackageSpec(url="https://github.com/lanl/JuliQAOA.jl")); Pkg.instantiate()'
julia --project=. -e 'import Pkg; Pkg.test(; coverage=false)'
```

The main CI workflow runs the package tests on Julia 1.10 and 1.11.

## Documentation

Build the Documenter site locally from the `docs/` environment:

```shell
julia --project=docs -e 'import Pkg; Pkg.develop([Pkg.PackageSpec(url="https://github.com/lanl/JuliQAOA.jl"), Pkg.PackageSpec(path=pwd())]); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated site is written to `docs/build/`.

## Angle-Transfer Fixture

QiskitOpt angle transfer is documented and tested with the saved fixture at
`test/fixtures/qaoa_regression.toml`. Do not add QiskitOpt as a package or test
dependency solely to validate angle transfer. Regenerate the fixture with:

```shell
julia --project=. test/fixtures/regenerate_qaoa_regression_fixture.jl
```

When changing the fixture, keep the reproducible QUBO case, normalization data,
exported parameters, and expected sample/value records together.

## Pull Requests

- Keep changes scoped to the issue or behavior under review.
- Add or update focused tests when changing solver behavior, metadata, or docs
  contracts.
- Do not weaken existing tests or CI checks to make a change pass.
- Document any local checks you could not run.
