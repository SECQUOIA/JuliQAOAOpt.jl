using Test
using Pkg
using TOML

using JuliQAOAOpt
using JuliQAOAOpt: MOI, QUBODrivers, QUBOTools

const SampleReads = QUBOTools.__moi_num_reads()

function julia_compat_spec()
    project = TOML.parsefile(joinpath(dirname(@__DIR__), "Project.toml"))
    return Pkg.Types.semver_spec(project["compat"]["julia"])
end

function build_two_variable_model(; final_reads = 25, number_of_reads = nothing, max_variables = 24)
    model = MOI.instantiate(JuliQAOAOpt.Optimizer; with_bridge_type = Float64)
    x, _ = MOI.add_constrained_variables(model, fill(MOI.ZeroOne(), 2))
    objective = MOI.ScalarQuadraticFunction{Float64}(
        MOI.ScalarQuadraticTerm{Float64}[
            MOI.ScalarQuadraticTerm{Float64}(2.0, x[1], x[2]),
        ],
        MOI.ScalarAffineTerm{Float64}[
            MOI.ScalarAffineTerm{Float64}(-1.0, x[1]),
            MOI.ScalarAffineTerm{Float64}(2.0, x[2]),
        ],
        0.0,
    )

    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(objective)}(), objective)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, JuliQAOAOpt.NumberOfLayers(), 1)
    MOI.set(model, JuliQAOAOpt.BasinHoppingIterations(), 1)
    MOI.set(model, JuliQAOAOpt.RandomSeed(), 11)
    MOI.set(model, JuliQAOAOpt.MaximumVariables(), max_variables)
    if !isnothing(number_of_reads)
        MOI.set(model, JuliQAOAOpt.NumberOfReads(), number_of_reads)
    end
    if !isnothing(final_reads)
        MOI.set(model, QUBODrivers.FinalNumberOfReads(), final_reads)
    end
    return model
end

function solution(model)
    return QUBOTools.solution(MOI.get(model, MOI.RawSolver()))
end

function qaoa_regression_fixture()
    return TOML.parsefile(joinpath(@__DIR__, "fixtures", "qaoa_regression.toml"))
end

@testset "support policy docs" begin
    readme = read(joinpath(dirname(@__DIR__), "README.md"), String)
    project = read(joinpath(dirname(@__DIR__), "Project.toml"), String)
    installation = findfirst("## Installation", readme)
    usage = findfirst("## Usage", readme)
    compat = julia_compat_spec()

    @test occursin("Julia 1.10 LTS and Julia 1.11", readme)
    @test occursin("JuliQAOA.find_angles_bh", readme)
    @test occursin("Julia 1.12+", readme)
    @test occursin("julia = \"1.10 - 1.11\"", project)
    @test VersionNumber("1.9.4") ∉ compat
    @test VersionNumber("1.10.0") ∈ compat
    @test VersionNumber("1.11.9") ∈ compat
    @test VersionNumber("1.12.0") ∉ compat
    @test occursin("Pkg.add(url=\"https://github.com/lanl/JuliQAOA.jl\")", readme)
    @test occursin("Pkg.add(url=\"https://github.com/SECQUOIA/JuliQAOAOpt.jl\")", readme)
    @test occursin("Pkg.add(\"JuMP\")", readme)
    @test installation !== nothing
    @test usage !== nothing

    if installation !== nothing && usage !== nothing
        @test first(installation) < first(usage)
    end
end

@testset "QUBODrivers generic tests" begin
    QUBODrivers.test(JuliQAOAOpt.Optimizer) do model
        MOI.set(model, MOI.Silent(), true)
        MOI.set(model, JuliQAOAOpt.NumberOfLayers(), 1)
        MOI.set(model, JuliQAOAOpt.BasinHoppingIterations(), 1)
        MOI.set(model, JuliQAOAOpt.RandomSeed(), 7)
        MOI.set(model, JuliQAOAOpt.MaximumVariables(), 8)
        MOI.set(model, QUBODrivers.FinalNumberOfReads(), 16)
    end
end

@testset "focused QAOA math regression fixture" begin
    fixture = qaoa_regression_fixture()
    energies = Float64.(fixture["qubo"]["energies"])
    angles = Float64.(fixture["parameters"]["normalized_angles"])
    layers = fixture["parameters"]["layers"]

    @test length(angles) == 2 * layers
    @test angles[1:layers] == [0.11, 0.22]
    @test angles[(layers + 1):end] == [0.33, 0.44]

    zscore = fixture["normalization"]["zscore"]
    normalized, shift, scale = JuliQAOAOpt._normalize_energies(energies, :zscore)

    @test shift ≈ zscore["shift"]
    @test scale ≈ zscore["scale"]
    @test normalized ≈ Float64.(zscore["normalized_energies"])
    @test JuliQAOAOpt._qiskit_parameters(angles, scale) ≈
          Float64.(fixture["parameters"]["zscore_qiskit_initial_parameters"])
    @test JuliQAOAOpt._qiskit_parameters(angles, scale)[1:layers] ≈ angles[1:layers]
    @test JuliQAOAOpt._qiskit_parameters(angles, scale)[(layers + 1):end] ≈
          angles[(layers + 1):end] ./ scale

    none = fixture["normalization"]["none"]
    unnormalized, none_shift, none_scale = JuliQAOAOpt._normalize_energies(energies, :none)

    @test none_shift ≈ none["shift"]
    @test none_scale ≈ none["scale"]
    @test unnormalized ≈ Float64.(none["normalized_energies"])
    @test JuliQAOAOpt._qiskit_parameters(angles, none_scale) ≈
          Float64.(fixture["parameters"]["none_qiskit_initial_parameters"])

    sampling = fixture["sampling"]
    draws = Float64.(sampling["draws"])
    samples = JuliQAOAOpt._sample_probabilities(
        Float64,
        Float64.(sampling["probabilities"]),
        energies,
        fixture["qubo"]["variables"],
        draws,
    )

    actual = Dict(
        Tuple(QUBOTools.state(sample)) =>
            (energy = QUBOTools.value(sample), reads = QUBOTools.reads(sample))
        for sample in samples
    )
    expected = Dict(
        Tuple(Int.(record["state"])) =>
            (energy = Float64(record["energy"]), reads = record["reads"])
        for record in sampling["samples"]
    )

    @test keys(actual) == keys(expected)
    for (state, record) in expected
        @test actual[state].energy ≈ record.energy
        @test actual[state].reads == record.reads
    end
    @test sum(QUBOTools.reads(sample) for sample in samples) == length(draws)
end

@testset "sampling and metadata" begin
    model = build_two_variable_model(final_reads = 31)
    MOI.optimize!(model)
    sampleset = solution(model)
    metadata = QUBOTools.metadata(sampleset)

    @test sum(QUBOTools.reads(sample) for sample in sampleset) == 31
    @test all(all(bit in (0, 1) for bit in QUBOTools.state(sample)) for sample in sampleset)

    n, L, Q, α, β = QUBOTools.qubo(MOI.get(model, MOI.RawSolver()), :dict; sense = :min)
    @test n == 2
    for sample in sampleset
        state = QUBOTools.state(sample)
        @test QUBOTools.value(sample) ≈ QUBOTools.value(state, L, Q, α, β)
    end

    @test haskey(metadata, "juliqaoa")
    data = metadata["juliqaoa"]
    @test data["number_of_layers"] == 1
    @test data["basinhopping_niter"] == 1
    @test data["enumerated_states"] == 4
    @test haskey(data, "energy_shift")
    @test haskey(data, "energy_scale")
    @test haskey(data, "expected_qubo_energy")
    @test length(data["normalized_angles"]) == 2
    @test length(data["qiskit_initial_parameters"]) == 2

    params = JuliQAOAOpt.qiskit_initial_parameters(sampleset)
    @test params == Float64.(data["qiskit_initial_parameters"])
    @test params[1] ≈ data["normalized_angles"][1]
    @test params[2] ≈ data["normalized_angles"][2] / data["energy_scale"]
end

@testset "metadata read counts and schema" begin
    expected_top_level_keys = Set([
        "algorithm",
        "backend",
        "execution",
        "juliqaoa",
        "optimizer",
        "origin",
        "reads",
        "seeds",
        "status",
        "termination_status",
        "time",
    ])

    default_model = build_two_variable_model(final_reads = nothing, number_of_reads = 7)
    MOI.optimize!(default_model)
    default_sampleset = solution(default_model)
    default_metadata = QUBOTools.metadata(default_sampleset)

    @test sum(QUBOTools.reads(sample) for sample in default_sampleset) == 7
    @test default_metadata["reads"]["number_of_reads"] == 7
    @test default_metadata["reads"]["final_number_of_reads"] == 7
    @test default_metadata["optimizer"]["iterations"] == 1
    @test default_metadata["optimizer"]["evaluations"] === nothing
    @test Set(keys(default_metadata)) == expected_top_level_keys
    @test Set(keys(default_metadata["time"])) == Set(["total"])

    default_data = default_metadata["juliqaoa"]
    @test default_data["enumerated_states"] == 4
    @test default_data["configured_number_of_reads"] == 7
    @test haskey(default_data, "time")
    @test !haskey(default_metadata, "enumerated_states")
    @test !haskey(default_metadata, "expected_qubo_energy")

    override_model = build_two_variable_model(final_reads = 5, number_of_reads = 13)
    MOI.optimize!(override_model)
    override_sampleset = solution(override_model)
    override_metadata = QUBOTools.metadata(override_sampleset)

    @test sum(QUBOTools.reads(sample) for sample in override_sampleset) == 5
    @test override_metadata["reads"]["number_of_reads"] == 5
    @test override_metadata["reads"]["final_number_of_reads"] == 5
    @test override_metadata["optimizer"]["evaluations"] === nothing
    @test Set(keys(override_metadata)) == expected_top_level_keys
    @test override_metadata["juliqaoa"]["configured_number_of_reads"] == 13
end

@testset "maximum variable guard" begin
    model = build_two_variable_model(max_variables = 1)
    @test_throws ArgumentError MOI.optimize!(model)
end
