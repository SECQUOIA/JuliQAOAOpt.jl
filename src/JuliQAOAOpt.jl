module JuliQAOAOpt

using JuliQAOA
using Random
using Statistics

import QUBODrivers:
    MOI,
    QUBODrivers,
    QUBOTools,
    Sample,
    SampleSet,
    sample

export Optimizer, qiskit_initial_parameters

QUBODrivers.@setup Optimizer begin
    name       = "JuliQAOA Optimizer"
    version    = v"0.1.0"
    attributes = begin
        NumberOfReads["num_reads"]::Integer                  = 1_000
        NumberOfLayers["num_layers"]::Integer                = 1
        BasinHoppingIterations["basinhopping_niter"]::Integer = 10
        RandomSeed["seed"]::Union{Integer,Nothing}           = nothing
        MaximumVariables["max_variables"]::Integer           = 24
        EnergyNormalization["energy_normalization"]::Symbol  = :zscore
    end
end

"""
    Optimizer

QUBODrivers-compatible MathOptInterface optimizer backed by JuliQAOA statevector
angle search and probability sampling.
"""
Optimizer

"""
    NumberOfReads()

Requested number of reads before QUBODrivers final-read normalization.
"""
NumberOfReads

"""
    NumberOfLayers()

Number of QAOA layers used for JuliQAOA angle search.
"""
NumberOfLayers

"""
    BasinHoppingIterations()

Number of basin-hopping iterations passed to `JuliQAOA.find_angles_bh`.
"""
BasinHoppingIterations

"""
    RandomSeed()

Optional non-negative integer seed used for angle search initialization and sampling.
"""
RandomSeed

"""
    MaximumVariables()

Maximum number of binary variables allowed before dense `2^n` state enumeration.
"""
MaximumVariables

"""
    EnergyNormalization()

Energy normalization policy used before angle search. Supported values are
`:zscore` and `:none`.
"""
EnergyNormalization

"""
    qiskit_initial_parameters(sampleset::SampleSet)
    qiskit_initial_parameters(metadata::AbstractDict)

Return QiskitOpt-ready QAOA initial parameters from JuliQAOAOpt sample-set metadata.
The order is beta values followed by gamma values, with gamma values converted from
normalized JuliQAOA energies back to the original QUBO energy scale.
"""
qiskit_initial_parameters(sampleset::SampleSet) = qiskit_initial_parameters(QUBOTools.metadata(sampleset))

function qiskit_initial_parameters(metadata::AbstractDict)
    data = get(metadata, "juliqaoa", nothing)
    data isa AbstractDict || throw(ArgumentError("metadata does not contain a 'juliqaoa' dictionary"))
    params = get(data, "qiskit_initial_parameters", nothing)
    isnothing(params) && throw(ArgumentError("metadata does not contain JuliQAOA qiskit_initial_parameters"))
    return Float64.(params)
end

function sample(sampler::Optimizer{T}) where {T}
    n, L, Q, α, β = QUBOTools.qubo(sampler, :dict; sense = :min)

    num_reads = MOI.get(sampler, NumberOfReads())
    final_reads = MOI.get(sampler, QUBODrivers.FinalNumberOfReads())
    p = MOI.get(sampler, NumberOfLayers())
    niter = MOI.get(sampler, BasinHoppingIterations())
    seed = MOI.get(sampler, RandomSeed())
    max_variables = MOI.get(sampler, MaximumVariables())
    normalization = MOI.get(sampler, EnergyNormalization())
    silent = MOI.get(sampler, MOI.Silent())

    _validate_inputs(n, num_reads, final_reads, p, niter, seed, max_variables, normalization)

    state_count = 1 << n
    enumeration = @timed _enumerate_energies(T, n, L, Q, α, β, state_count)
    energies = enumeration.value
    normalized, shift, scale = _normalize_energies(energies, normalization)

    if !isnothing(seed)
        Random.seed!(seed)
    end

    mixer = JuliQAOA.mixer_x(n)
    angle_search = @timed JuliQAOA.find_angles_bh(
        p,
        mixer,
        normalized;
        max = false,
        niter = niter,
        verbose = !silent,
    )
    angle_sets, expected_normalized_values = angle_search.value
    normalized_angles = Float64.(angle_sets[p])

    probability_result = @timed _probabilities(normalized_angles, mixer, normalized)
    probabilities = probability_result.value
    expected_energy = float(sum(probabilities .* Float64.(energies)))
    expected_normalized_energy = float(sum(probabilities .* normalized))

    sampling = @timed _sample_probabilities(T, probabilities, energies, n, final_reads, seed)
    samples = sampling.value
    effective_time = enumeration.time + angle_search.time + probability_result.time + sampling.time

    metadata = QUBODrivers._sampler_metadata(
        origin                = "JuliQAOAOpt.jl using JuliQAOA.jl",
        algorithm_name        = "JuliQAOA",
        execution_mode        = "statevector_angle_search",
        optimizer_iterations  = niter,
        number_of_reads       = final_reads,
        final_number_of_reads = final_reads,
        seeds                 = Dict{String,Any}("sampler" => seed),
        status                = "locally_solved",
        termination_status    = MOI.LOCALLY_SOLVED,
    )
    metadata["time"] = Dict{String,Any}("effective" => effective_time)
    metadata["juliqaoa"] = Dict{String,Any}(
        "number_of_layers" => p,
        "basinhopping_niter" => niter,
        "configured_number_of_reads" => num_reads,
        "energy_normalization" => string(normalization),
        "energy_shift" => shift,
        "energy_scale" => scale,
        "enumerated_states" => state_count,
        "expected_qubo_energy" => expected_energy,
        "expected_normalized_energy" => expected_normalized_energy,
        "expected_values_by_layer" => Float64.(expected_normalized_values),
        "normalized_angles" => normalized_angles,
        "qiskit_initial_parameters" => _qiskit_parameters(normalized_angles, scale),
        "time" => Dict{String,Any}(
            "enumeration" => enumeration.time,
            "angle_search" => angle_search.time,
            "probabilities" => probability_result.time,
            "sampling" => sampling.time,
            "effective" => effective_time,
        ),
    )

    return SampleSet{T}(samples, metadata; sense = :min, domain = :bool)
end

function _validate_inputs(n, num_reads, final_reads, p, niter, seed, max_variables, normalization)
    n >= 1 || throw(ArgumentError("JuliQAOAOpt requires at least one variable"))
    n <= max_variables ||
        throw(ArgumentError("QUBO has $(n) variables, exceeding MaximumVariables()=$(max_variables)"))
    num_reads >= 0 || throw(ArgumentError("num_reads must be non-negative"))
    final_reads >= 0 || throw(ArgumentError("final_num_reads must be non-negative"))
    p >= 1 || throw(ArgumentError("num_layers must be at least 1"))
    niter >= 0 || throw(ArgumentError("basinhopping_niter must be non-negative"))
    isnothing(seed) || seed >= 0 || throw(ArgumentError("seed must be non-negative or nothing"))
    normalization in (:zscore, :none) || throw(ArgumentError("energy_normalization must be :zscore or :none"))
    return nothing
end

function _enumerate_energies(::Type{T}, n, L, Q, α, β, state_count) where {T}
    energies = Vector{T}(undef, state_count)
    for i in 1:state_count
        state = _state(i, n)
        energies[i] = QUBOTools.value(state, L, Q, α, β)
    end
    return energies
end

_state(index::Integer, n::Integer) = digits(Int, index - 1; base = 2, pad = n)

function _normalize_energies(energies, normalization::Symbol)
    values = Float64.(energies)
    if normalization == :none
        return values, 0.0, 1.0
    end

    shift = mean(values)
    scale = std(values)
    if !isfinite(scale) || iszero(scale)
        scale = 1.0
    end
    return (values .- shift) ./ scale, shift, scale
end

function _probabilities(angles, mixer, normalized)
    values = Float64.(real.(JuliQAOA.probabilities(angles, mixer, normalized)))
    values .= max.(values, 0.0)
    total = sum(values)
    if total <= 0 || !isfinite(total)
        throw(ErrorException("JuliQAOA returned invalid probabilities"))
    end
    return values ./ total
end

function _sample_probabilities(::Type{T}, probabilities, energies, n, final_reads, seed) where {T}
    rng = isnothing(seed) ? Random.default_rng() : Random.MersenneTwister(seed)
    draws = rand(rng, final_reads)
    return _sample_probabilities(T, probabilities, energies, n, draws)
end

function _sample_probabilities(
    ::Type{T},
    probabilities,
    energies,
    n,
    draws::AbstractVector{<:Real},
) where {T}
    cumulative = cumsum(probabilities)
    cumulative[end] = 1.0
    counts = Dict{Int,Int}()
    for u in draws
        index = first(searchsorted(cumulative, u))
        counts[index] = get(counts, index, 0) + 1
    end

    samples = Sample{T,Int}[]
    sizehint!(samples, length(counts))
    for (index, reads) in counts
        push!(samples, Sample{T,Int}(_state(index, n), energies[index], reads))
    end
    return samples
end

function _qiskit_parameters(angles, scale)
    p = length(angles) ÷ 2
    return Float64[v for v in vcat(angles[1:p], angles[(p + 1):end] ./ scale)]
end

end # module
