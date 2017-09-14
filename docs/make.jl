using Documenter, MIToS

# TO DO: ADD StatPlots to DOC REQUIRES

makedocs(
    format = :html,
    sitename = "MIToS",
    modules = [MIToS],
    pages = [
        "man/Installation.md",
        "man/Example.md",
        "Modules" => [  "man/MSA.md",
                        "man/Information.md",
                        "man/SIFTS.md",
                        "man/PDB.md",
                        "man/Pfam.md"    ],
        "man/Scripts.md",
        "API" => [  "man/MSA_API.md",
                    "man/Information_API.md",
                    "man/SIFTS_API.md",
                    "man/PDB_API.md",
                    "man/Pfam_API.md",
                    "man/Utils_API.md"          ]
    ]
)

deploydocs(
    repo   = "github.com/diegozea/MIToS.jl.git",
    target = "build",
    julia  = "0.5",
    deps   = nothing,
    make   = nothing
)
