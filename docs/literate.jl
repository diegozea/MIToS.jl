using Literate, MIToS

COOKBOOK = joinpath(@__DIR__, "src", "cookbook")
MD_OUTPUT = joinpath(@__DIR__, "src")
NB_OUTPUT = joinpath(@__DIR__, "src", "cookbook", "notebooks")

for file in [ "01_Change_B_factors.jl",
              "02_Linking_structural_and_evolutionary_information.jl"]
    Literate.markdown(joinpath(COOKBOOK, file),
        MD_OUTPUT, execute=false, documenter=true)
    Literate.notebook(joinpath(COOKBOOK, file),
        NB_OUTPUT, execute=false, documenter=true)
end
