using Documenter, MIToS

makedocs(
    doctest = true,
    format = :html,
    sitename = "MIToS",
    modules = [MIToS],
    pages = [
        "Home" => "index.md",
        "Installation.md",
        "Example.md",
        "Modules" => [  "MSA.md",
                        "Information.md",
                        "SIFTS.md",
                        "PDB.md",
                        "Pfam.md"    ],
        "Scripts.md",
        "API" => [  "MSA_API.md",
                    "Information_API.md",
                    "SIFTS_API.md",
                    "PDB_API.md",
                    "Pfam_API.md",
                    "Utils_API.md"          ]
    ]
)

deploydocs(
    repo   = "github.com/diegozea/MIToS.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
