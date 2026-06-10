using Documenter
using JuliQAOAOpt

DocMeta.setdocmeta!(JuliQAOAOpt, :DocTestSetup, :(using JuliQAOAOpt); recursive = true)

makedocs(
    modules = [JuliQAOAOpt],
    sitename = "JuliQAOAOpt.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Usage" => "manual/usage.md",
            "Metadata" => "manual/metadata.md",
            "QiskitOpt angle transfer" => "manual/qiskit.md",
        ],
        "API" => "api.md",
    ],
    checkdocs = :none,
)
