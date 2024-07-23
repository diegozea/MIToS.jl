using Documenter, MIToS

const WARNONLY = [:missing_docs]
if get(ENV, "CI", nothing) === nothing
    @info "Running locally, adding :cross_references to WARNONLY"
    push!(WARNONLY, :cross_references)
end

DocMeta.setdocmeta!(MIToS, :DocTestSetup, :(using MIToS); recursive = true)

include("literate.jl")

makedocs(
    doctest = true,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "MIToS",
    authors = "Diego Javier Zea",
    modules = [MIToS],
    pages = [
        "Home" => "index.md",
        "Installation.md",
        "Example.md",
        "Modules" => ["MSA.md", "Information.md", "SIFTS.md", "PDB.md", "Pfam.md"],
        "Cookbook" => [
            "01_Change_B_factors.md",
            "02_Linking_structural_and_evolutionary_information.md",
            "03_RMSF.md",
        ],
        "API" => [
            "MSA_API.md",
            "Information_API.md",
            "SIFTS_API.md",
            "PDB_API.md",
            "Pfam_API.md",
            "Utils_API.md",
        ],
        "Scripts.md",
    ],
    warnonly = WARNONLY,
)

deploydocs(
    repo = "github.com/diegozea/MIToS.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
