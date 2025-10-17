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
project_folder = normpath(@__DIR__, "..")
Pkg.activate(project_folder)
Pkg.instantiate()
using MIToS_Scripts
using MIToS
scripts_folder = joinpath(pkgdir(MIToS_Scripts), "scripts")

function run_script_help(script_name)
    script_path = joinpath(scripts_folder, script_name)
    cmd = `$(Base.julia_cmd()) --project=$project_folder -e 'using MIToS_Scripts, MIToS; MIToS_Scripts.loadedversion(::Module) = Base.pkgversion(MIToS); script_path = popfirst!(Base.ARGS); include(script_path)' $script_path -h`
    read(cmd, String)
end
```

### Buslje09.jl

```@example scripts
script_name = "Buslje09.jl" # hide
print(run_script_help(script_name)) # hide
```

### BLMI.jl

```@example scripts
script_name = "BLMI.jl" # hide
print(run_script_help(script_name)) # hide

```

### Conservation.jl

```@example scripts
script_name = "Conservation.jl" # hide
print(run_script_help(script_name)) # hide
```

### DownloadPDB.jl

```@example scripts
script_name = "DownloadPDB.jl" # hide
print(run_script_help(script_name)) # hide
```

### Distances.jl

```@example scripts
script_name = "Distances.jl" # hide
print(run_script_help(script_name)) # hide
```

### MSADescription.jl

```@example scripts
script_name = "MSADescription.jl" # hide
print(run_script_help(script_name)) # hide
```

### PercentIdentity.jl

```@example scripts
script_name = "PercentIdentity.jl" # hide
print(run_script_help(script_name)) # hide
```

### AlignedColumns.jl

```@example scripts
script_name = "AlignedColumns.jl" # hide
print(run_script_help(script_name)) # hide
```

### SplitStockholm.jl

```@example scripts
script_name = "SplitStockholm.jl" # hide
print(run_script_help(script_name)) # hide
```
