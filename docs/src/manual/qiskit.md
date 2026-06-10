# QiskitOpt Angle Transfer

JuliQAOAOpt stores learned normalized angles and the converted QiskitOpt initial
parameter vector in sample-set metadata. The QiskitOpt vector is beta values first,
then gamma values. Gamma values are rescaled from the normalized JuliQAOA energy
diagonal back to the original QUBO energy scale.

```julia
using JuliQAOAOpt
using QUBOTools

sampleset = QUBOTools.solution(unsafe_backend(model))
qiskit_params = JuliQAOAOpt.qiskit_initial_parameters(sampleset)
```

Pass that vector to QiskitOpt when building the same QUBO there:

```julia
using QiskitOpt
using QiskitOpt: QAOA

qiskit_model = Model(QiskitOpt.QAOA.Optimizer)
# Build the same QUBO on qiskit_model.

set_attribute(qiskit_model, QAOA.NumberOfLayers(), 2)
set_attribute(qiskit_model, QAOA.InitialParameters(), qiskit_params)
set_attribute(qiskit_model, QAOA.InitialParameterSource(), "JuliQAOAOpt")
```

QiskitOpt is intentionally not a JuliQAOAOpt package or test dependency. The
repository keeps a small saved fixture at `test/fixtures/qaoa_regression.toml`
instead. That fixture records:

- a reproducible two-variable QUBO, `-x1 + 2*x2 + 2*x1*x2`;
- all Boolean states and expected QUBO values;
- `:zscore` and `:none` normalization data;
- exported normalized angles and QiskitOpt initial parameters;
- deterministic sample records generated from saved probability draws.

For the fixture's two-layer angle vector `[0.11, 0.22, 0.33, 0.44]`, the saved
z-score transfer data expects:

```julia
[0.11, 0.22, 0.1807484439767048, 0.24099792530227307]
```

The fixture can be regenerated with:

```shell
julia --project=. test/fixtures/regenerate_qaoa_regression_fixture.jl
```
