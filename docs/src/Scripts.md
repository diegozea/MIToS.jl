```@setup log
@info "Scripts docs"
```

# Scripts

MIToS implements several useful scripts to **command line execution
(without requiring Julia coding)**. All this scripts are located in the `scripts` folder
of the MIToS directory. You can copy them to your working directory, use the path to
their folder or put them in the path
(look into the **Installation** section of this manual).  

```@contents
Pages = ["Scripts.md"]
Depth = 4
```   

```@setup
julia = Base.julia_cmd(); # path to the julia executable
run(`$julia -e 'using Pkg; Pkg.add("MIToS")'`) # add MIToS
```

## Buslje09.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
script_path = joinpath(pkgdir(MIToS), "scripts", "Buslje09.jl")
run(`$julia $script_path -h`)
```  

## BLMI.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
script_path = joinpath(pkgdir(MIToS), "scripts", "BLMI.jl")
run(`$julia $script_path -h`)
```  

## Conservation.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
script_path = joinpath(pkgdir(MIToS), "scripts", "Conservation.jl")
run(`$julia $script_path -h`)
```  

## DownloadPDB.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
script_path = joinpath(pkgdir(MIToS), "scripts", "DownloadPDB.jl")
run(`$julia $script_path -h`)
```  

## Distances.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
script_path = joinpath(pkgdir(MIToS), "scripts", "Distances.jl")
run(`$julia $script_path -h`)
```  

## MSADescription.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
script_path = joinpath(pkgdir(MIToS), "scripts", "MSADescription.jl")
run(`$julia $script_path -h`)
```  

## PercentIdentity.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
script_path = joinpath(pkgdir(MIToS), "scripts", "PercentIdentity.jl")
run(`$julia $script_path -h`)
```  

## AlignedColumns.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
script_path = joinpath(pkgdir(MIToS), "scripts", "AlignedColumns.jl")
run(`$julia $script_path -h`)
```  

## SplitStockholm.jl

```@repl
using MIToS
julia = Base.julia_cmd(); # path to the julia executable
script_path = joinpath(pkgdir(MIToS), "scripts", "SplitStockholm.jl")
run(`$julia $script_path -h`)
```

