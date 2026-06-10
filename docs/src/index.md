# JuliQAOAOpt.jl

`JuliQAOAOpt.jl` is a JuMP and QUBODrivers interface for LANL's
[`JuliQAOA.jl`](https://github.com/lanl/JuliQAOA.jl). It converts QUBO models into
a local statevector QAOA problem, learns angles with `JuliQAOA.find_angles_bh`,
samples from the resulting exact probabilities, and returns a QUBODrivers-compatible
sample set.

This package is intended for moderate-size QUBOs. It enumerates all `2^n` Boolean
states before angle search, so the `MaximumVariables()` attribute is an important
local guardrail.

## Package Scope

These docs describe the JuliQAOA-specific interface: installation, optimizer
attributes, metadata, scalability limits, and QiskitOpt angle transfer.

For the generic mechanics of implementing a QUBODrivers sampler, use the
[QUBODrivers sampler setup guide](https://juliaqubo.github.io/QUBODrivers.jl/stable/manual/4-setup/)
instead of duplicating that tutorial here.

## Contents

```@contents
Pages = [
    "manual/usage.md",
    "manual/metadata.md",
    "manual/qiskit.md",
    "api.md",
]
Depth = 2
```
