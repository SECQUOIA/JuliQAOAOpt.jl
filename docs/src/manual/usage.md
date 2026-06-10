# Usage

## Installation

JuliQAOAOpt supports Julia 1.10 LTS and Julia 1.11 for angle finding. The optimizer
calls `JuliQAOA.find_angles_bh`, which currently relies on JuliQAOA and Enzyme
gradient computation that upstream JuliQAOA documents as limited on Julia 1.12+.
Until that path is reliable on Julia 1.12+, this package constrains Julia compat to
`1.10 - 1.11`.

`JuliQAOA.jl` is not currently registered in Julia's General registry, so install it
by URL before installing this package:

```julia
import Pkg
Pkg.add(url = "https://github.com/lanl/JuliQAOA.jl")
Pkg.add(url = "https://github.com/SECQUOIA/JuliQAOAOpt.jl")
Pkg.add("JuMP")
```

## Solving A QUBO

```julia
using JuMP
using JuliQAOAOpt

model = Model(JuliQAOAOpt.Optimizer)

Q = [-1  2  2
      2 -1  2
      2  2 -1]

@variable(model, x[1:3], Bin)
@objective(model, Min, x' * Q * x)

set_attribute(model, JuliQAOAOpt.NumberOfLayers(), 2)
set_attribute(model, JuliQAOAOpt.BasinHoppingIterations(), 10)
set_attribute(model, JuliQAOAOpt.RandomSeed(), 1234)
set_attribute(model, JuliQAOAOpt.MaximumVariables(), 24)

optimize!(model)
```

## Optimizer Attributes

| Attribute | Storage key | Default | Meaning |
| --- | --- | ---: | --- |
| `NumberOfReads()` | `num_reads` | `1000` | Requested reads before QUBODrivers final-read normalization. |
| `NumberOfLayers()` | `num_layers` | `1` | QAOA layer count passed to JuliQAOA angle search. |
| `BasinHoppingIterations()` | `basinhopping_niter` | `10` | Basin-hopping iterations for `JuliQAOA.find_angles_bh`. |
| `RandomSeed()` | `seed` | `nothing` | Optional non-negative seed for angle search initialization and sampling. |
| `MaximumVariables()` | `max_variables` | `24` | Guardrail checked before dense `2^n` state enumeration. |
| `EnergyNormalization()` | `energy_normalization` | `:zscore` | Energy normalization before angle search. Supported values are `:zscore` and `:none`. |

`QUBODrivers.FinalNumberOfReads()` can override the final number of returned sample
reads. When it is not set, the configured `NumberOfReads()` value is used.

## Scalability Limits

JuliQAOAOpt builds JuliQAOA's cost diagonal by enumerating every Boolean state for
an `n`-variable QUBO. Runtime and memory therefore grow exponentially with `n`.

The default `MaximumVariables()` value is `24`, corresponding to
`2^24 = 16,777,216` states. A single dense `Float64` vector at that size is about
128 MiB, and the full solve allocates additional normalized energy, probability,
and statevector data.

If a model exceeds the configured cap, `optimize!` throws an `ArgumentError` before
enumerating states:

```julia
set_attribute(model, JuliQAOAOpt.MaximumVariables(), 20)
```
