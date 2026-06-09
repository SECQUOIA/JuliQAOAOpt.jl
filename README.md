# JuliQAOAOpt.jl
[![QUBODRIVERS](https://img.shields.io/badge/Powered%20by-QUBODrivers.jl-%20%234063d8)](https://github.com/JuliaQUBO/QUBODrivers.jl)

JuMP/QUBODrivers interface for LANL's [JuliQAOA.jl](https://github.com/lanl/JuliQAOA.jl).

`JuliQAOAOpt.jl` converts QUBO models into a local statevector QAOA problem, learns QAOA
angles with `JuliQAOA.find_angles_bh`, samples from the resulting exact probabilities, and
returns a QUBODrivers-compatible sample set.

This is intended for moderate-size QUBOs. The current implementation enumerates all
``2^n`` Boolean states to build the JuliQAOA cost diagonal.

## Installation

JuliQAOAOpt supports Julia 1.10 LTS and Julia 1.11 for angle finding. The optimizer calls
`JuliQAOA.find_angles_bh`, which currently relies on JuliQAOA/Enzyme gradient computation
that upstream JuliQAOA documents as limited on Julia 1.12+. Until that path is reliable on
Julia 1.12+, this package constrains Julia compat to `1.10 - 1.11`.

`JuliQAOA.jl` is not currently registered in Julia's General registry, so install it by URL
before installing this package:

```julia
import Pkg
Pkg.add(url="https://github.com/lanl/JuliQAOA.jl")
Pkg.add(url="https://github.com/SECQUOIA/JuliQAOAOpt.jl")
Pkg.add("JuMP")
```

## Usage

```julia
using JuMP
using JuliQAOAOpt

model = Model(JuliQAOAOpt.Optimizer)

Q = [ -1  2  2
       2 -1  2
       2  2 -1 ]

@variable(model, x[1:3], Bin)
@objective(model, Min, x' * Q * x)

set_attribute(model, JuliQAOAOpt.NumberOfLayers(), 2)
set_attribute(model, JuliQAOAOpt.BasinHoppingIterations(), 10)
set_attribute(model, JuliQAOAOpt.RandomSeed(), 1234)
set_attribute(model, JuliQAOAOpt.MaximumVariables(), 24)

optimize!(model)
```

`JuliQAOAOpt` stores the learned normalized angles, energy normalization, and QiskitOpt-ready
parameters in the returned sample-set metadata.

## Transferring Angles To QiskitOpt.jl

```julia
using QiskitOpt
using QiskitOpt: QAOA
using QUBOTools

sampleset = QUBOTools.solution(unsafe_backend(model))
qiskit_params = JuliQAOAOpt.qiskit_initial_parameters(sampleset)

qiskit_model = Model(QiskitOpt.QAOA.Optimizer)
# build the same QUBO on qiskit_model

set_attribute(qiskit_model, QAOA.NumberOfLayers(), 2)
set_attribute(qiskit_model, QAOA.InitialParameters(), qiskit_params)
set_attribute(qiskit_model, QAOA.InitialParameterSource(), "JuliQAOAOpt")
```

The returned parameters are in Qiskit's beta-then-gamma order. Gamma values are rescaled from
the normalized JuliQAOA energy diagonal back to the original QUBO energy scale.
