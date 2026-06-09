using Statistics

function state(index::Integer, n::Integer)
    return digits(Int, index - 1; base = 2, pad = n)
end

format_float(value) = repr(Float64(value))
format_int_array(values) = "[" * join(string.(values), ", ") * "]"
format_float_array(values) = "[" * join(format_float.(values), ", ") * "]"
format_state_array(states) = "[" * join(format_int_array.(states), ", ") * "]"

function qiskit_parameters(angles, scale)
    layers = length(angles) ÷ 2
    return Float64[v for v in vcat(angles[1:layers], angles[(layers + 1):end] ./ scale)]
end

function sample_counts(probabilities, draws)
    cumulative = cumsum(probabilities)
    cumulative[end] = 1.0
    counts = Dict{Int,Int}()
    for u in draws
        index = first(searchsorted(cumulative, u))
        counts[index] = get(counts, index, 0) + 1
    end
    return sort(collect(counts); by = first)
end

fixture_path = joinpath(@__DIR__, "qaoa_regression.toml")

n = 2
states = [state(index, n) for index in 1:(1 << n)]
energies = [0.0, -1.0, 2.0, 3.0]
shift = mean(energies)
scale = std(energies)
normalized = (energies .- shift) ./ scale
angles = [0.11, 0.22, 0.33, 0.44]
probabilities = [0.05, 0.15, 0.30, 0.50]
draws = [0.01, 0.04, 0.06, 0.10, 0.15, 0.19, 0.21, 0.25, 0.30, 0.40, 0.45, 0.49, 0.51, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99, 0.999]
counts = sample_counts(probabilities, draws)

open(fixture_path, "w") do io
    println(io, "[fixture]")
    println(io, "description = \"Two-variable QUBO used for focused JuliQAOAOpt regression tests.\"")
    println(io, "regenerate = \"julia --project=. test/fixtures/regenerate_qaoa_regression_fixture.jl\"")
    println(io)
    println(io, "[qubo]")
    println(io, "variables = $n")
    println(io, "objective = \"-x1 + 2*x2 + 2*x1*x2\"")
    println(io, "states = $(format_state_array(states))")
    println(io, "energies = $(format_float_array(energies))")
    println(io)
    println(io, "[normalization.zscore]")
    println(io, "shift = $(format_float(shift))")
    println(io, "scale = $(format_float(scale))")
    println(io, "normalized_energies = $(format_float_array(normalized))")
    println(io)
    println(io, "[normalization.none]")
    println(io, "shift = 0.0")
    println(io, "scale = 1.0")
    println(io, "normalized_energies = $(format_float_array(energies))")
    println(io)
    println(io, "[parameters]")
    println(io, "layers = $(length(angles) ÷ 2)")
    println(io, "normalized_angles = $(format_float_array(angles))")
    println(io, "zscore_qiskit_initial_parameters = $(format_float_array(qiskit_parameters(angles, scale)))")
    println(io, "none_qiskit_initial_parameters = $(format_float_array(qiskit_parameters(angles, 1.0)))")
    println(io)
    println(io, "[sampling]")
    println(io, "draws = $(format_float_array(draws))")
    println(io, "probabilities = $(format_float_array(probabilities))")

    for (index, reads) in counts
        println(io)
        println(io, "[[sampling.samples]]")
        println(io, "state = $(format_int_array(states[index]))")
        println(io, "energy = $(format_float(energies[index]))")
        println(io, "reads = $reads")
    end
end
