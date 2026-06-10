using Documenter
using JuliQAOAOpt

DocMeta.setdocmeta!(JuliQAOAOpt, :DocTestSetup, :(using JuliQAOAOpt); recursive = true)

makedocs(
    modules = [JuliQAOAOpt],
    sitename = "JuliQAOAOpt.jl",
    clean = true,
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
    workdir = @__DIR__,
    checkdocs = :none,
)

if "--deploy" in ARGS
    deploydocs(repo = "github.com/SECQUOIA/JuliQAOAOpt.jl.git", push_preview = true)
else
    @warn "Skipping deployment"
end
