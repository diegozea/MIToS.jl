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

```@example scripts
using Pkg
Pkg.add(url="https://github.com/MIToSOrg/MIToS_Scripts.jl")
```

Then, you can get the location of the installed scripts by running the following command:

```@example scripts
using MIToS_Scripts
scripts_folder = joinpath(pkgdir(MIToS_Scripts), "scripts")
```

You can run them from that location or copy them to a directory in your `PATH`.

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

### Buslje09.jl

```@example scripts
script_path = joinpath(scripts_folder, "Buslje09.jl") # path to the script
run(`$(Base.julia_cmd()) $script_path -h`)
```  

### BLMI.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
scripts_folder = joinpath(pkgdir(MIToS), "scripts")
script_path = joinpath(scripts_folder, "BLMI.jl")
run(`$julia --project=$scripts_folder $script_path -h`)
```  

### Conservation.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
scripts_folder = joinpath(pkgdir(MIToS), "scripts")
script_path = joinpath(scripts_folder, "Conservation.jl")
run(`$julia --project=$scripts_folder $script_path -h`)
```  

### DownloadPDB.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
scripts_folder = joinpath(pkgdir(MIToS), "scripts")
script_path = joinpath(scripts_folder, "DownloadPDB.jl")
run(`$julia --project=$scripts_folder $script_path -h`)
```  

### Distances.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
scripts_folder = joinpath(pkgdir(MIToS), "scripts")
script_path = joinpath(scripts_folder, "Distances.jl")
run(`$julia --project=$scripts_folder $script_path -h`)
```  

### MSADescription.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
scripts_folder = joinpath(pkgdir(MIToS), "scripts")
script_path = joinpath(scripts_folder, "MSADescription.jl")
run(`$julia --project=$scripts_folder $script_path -h`)
```  

### PercentIdentity.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
scripts_folder = joinpath(pkgdir(MIToS), "scripts")
script_path = joinpath(scripts_folder, "PercentIdentity.jl")
run(`$julia --project=$scripts_folder $script_path -h`)
```  

### AlignedColumns.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
scripts_folder = joinpath(pkgdir(MIToS), "scripts")
script_path = joinpath(scripts_folder, "AlignedColumns.jl")
run(`$julia --project=$scripts_folder $script_path -h`)
```  

### SplitStockholm.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
scripts_folder = joinpath(pkgdir(MIToS), "scripts")
script_path = joinpath(scripts_folder, "SplitStockholm.jl")
run(`$julia --project=$scripts_folder $script_path -h`)
```

