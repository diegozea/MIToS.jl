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
                        "man/PDB.md"    ]
    ]
)

deploydocs(
    repo   = "github.com/diegozea/MIToS.jl.git",
    target = "build",
    julia  = "0.5",
    deps   = nothing,
    make   = nothing
)
