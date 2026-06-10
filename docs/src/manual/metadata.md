# Metadata

JuliQAOAOpt returns a QUBODrivers `SampleSet`. Standard QUBODrivers metadata remains
at the top level, while solver-specific details are stored under
`metadata["juliqaoa"]`.

```julia
using QUBOTools

sampleset = QUBOTools.solution(unsafe_backend(model))
metadata = QUBOTools.metadata(sampleset)
juliqaoa = metadata["juliqaoa"]
```

## `metadata["juliqaoa"]`

| Key | Meaning |
| --- | --- |
| `number_of_layers` | QAOA layer count used for the solve. |
| `basinhopping_niter` | Basin-hopping iterations requested for angle search. |
| `configured_number_of_reads` | Value from `NumberOfReads()` before final-read normalization. |
| `energy_normalization` | Normalization policy, currently `"zscore"` or `"none"`. |
| `energy_shift` | Additive shift used for normalized energies. |
| `energy_scale` | Multiplicative scale used for normalized energies. |
| `enumerated_states` | Number of Boolean states enumerated, equal to `2^n`. |
| `expected_qubo_energy` | Probability-weighted expected value on the original QUBO scale. |
| `expected_normalized_energy` | Probability-weighted expected value on the normalized energy scale. |
| `expected_values_by_layer` | Expected normalized values returned by JuliQAOA for each layer. |
| `normalized_angles` | Learned JuliQAOA angles on the normalized energy scale. |
| `qiskit_initial_parameters` | QiskitOpt-ready beta-then-gamma parameter vector. |
| `time` | Timing details for enumeration, angle search, probabilities, and sampling. |

Use [`JuliQAOAOpt.qiskit_initial_parameters`](@ref) to read the parameter vector
instead of reaching into the dictionary directly.
