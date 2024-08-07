```@setup log
@info "Scripts docs"
```

# MIToS' Scripts

The [MIToS_Scripts.jl](https://github.com/MIToSOrg/MIToS_Scripts.jl) package offers a set
of easy-to-use scripts for command-line execution without requiring Julia coding.
It includes several scripts designed for various bioinformatics tasks, such as measuring
estimating residue conservation and inter-residue coevolution, calculating distances between
residues in a protein structure, and more.

```@contents
Pages = ["Scripts.md"]
Depth = 4
```

## Installation

To install **MIToS_Scripts.jl**, you only need Julia 1.9 or later installed on your
system. Executing `julia` in the terminal to open the Julia REPL, and finally, run the
following command:

```julia
using Pkg
Pkg.add(url = "https://github.com/MIToSOrg/MIToS_Scripts.jl")
```

Then, you can get the location of the installed scripts by running the following command:

```julia
using MIToS_Scripts
scripts_folder = joinpath(pkgdir(MIToS_Scripts), "scripts")
```

You can run them from that location. Alternatively, you can add the location to your
`PATH` environment variable, or copy the scripts to a folder already in your `PATH` to
run them from anywhere.

## Usage

You can execute each provided script from your command line. For example, to run the `Buslje09.jl`
script—if you are located in the folder where it is the scripts—use:

```bash
julia Buslje09.jl input_msa_file
```

Refer to the documentation of each script for specific usage instructions; you can access
it by running the script with the `--help` or `-h` flag:

```bash
julia Buslje09.jl -h
```

## Scripts

```@setup scripts
using Pkg
project_folder = "MIToS_Scripts_Project"
isdir(project_folder) || mkdir(project_folder)
Pkg.activate(project_folder)
Pkg.add(url="https://github.com/MIToSOrg/MIToS_Scripts.jl")
using MIToS_Scripts
scripts_folder = joinpath(pkgdir(MIToS_Scripts), "scripts")
```

### Buslje09.jl

```@example scripts
script_path = joinpath(scripts_folder, "Buslje09.jl") # hide
println(read(`$(Base.julia_cmd()) --project=$project_folder $script_path -h`, String)) #hide
```

### BLMI.jl

```@example scripts
script_path = joinpath(scripts_folder, "BLMI.jl") # hide
println(read(`$(Base.julia_cmd()) --project=$project_folder $script_path -h`, String)) #hide

```

### Conservation.jl

```@example scripts
script_path = joinpath(scripts_folder, "Conservation.jl") # hide
println(read(`$(Base.julia_cmd()) --project=$project_folder $script_path -h`, String)) #hide
```

### DownloadPDB.jl

```@example scripts
script_path = joinpath(scripts_folder, "DownloadPDB.jl") # hide
println(read(`$(Base.julia_cmd()) --project=$project_folder $script_path -h`, String)) #hide
```

### Distances.jl

```@example scripts
script_path = joinpath(scripts_folder, "Distances.jl") # hide
println(read(`$(Base.julia_cmd()) --project=$project_folder $script_path -h`, String)) #hide
```

### MSADescription.jl

```@example scripts
script_path = joinpath(scripts_folder, "MSADescription.jl") # hide
println(read(`$(Base.julia_cmd()) --project=$project_folder $script_path -h`, String)) #hide
```

### PercentIdentity.jl

```@example scripts
script_path = joinpath(scripts_folder, "PercentIdentity.jl") # hide
println(read(`$(Base.julia_cmd()) --project=$project_folder $script_path -h`, String)) #hide
```

### AlignedColumns.jl

```@example scripts
script_path = joinpath(scripts_folder, "AlignedColumns.jl") # hide
println(read(`$(Base.julia_cmd()) --project=$project_folder $script_path -h`, String)) #hide
```

### SplitStockholm.jl

```@example scripts
script_path = joinpath(scripts_folder, "SplitStockholm.jl") # hide
println(read(`$(Base.julia_cmd()) --project=$project_folder $script_path -h`, String)) #hide
```
