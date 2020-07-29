using Documenter, MIToS

include("literate.jl")

makedocs(
    doctest = true,
	format = Documenter.HTML(),
    sitename = "MIToS",
    authors = "Diego Javier Zea",
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
        "Cookbook" => [ "01_Change_B_factors.md",
                        "02_Linking_structural_and_evolutionary_information.md",
						"03_RMSF.md"],
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
