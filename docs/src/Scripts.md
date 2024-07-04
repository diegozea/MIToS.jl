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

## Buslje09.jl

```@repl
using MIToS
julia = Base.julia_cmd()
script_path = joinpath(pkgdir(MIToS), "scripts", "Buslje09.jl")
script_output = read(`$julia $script_path -h`, String)
print(script_output)
```  

## BLMI.jl

```@repl
using MIToS
julia = Base.julia_cmd()
script_path = joinpath(pkgdir(MIToS), "scripts", "BLMI.jl")
script_output = read(`$julia $script_path -h`, String)
print(script_output)
```  

## Conservation.jl

```@repl
using MIToS
julia = Base.julia_cmd()
script_path = joinpath(pkgdir(MIToS), "scripts", "Conservation.jl")
script_output = read(`$julia $script_path -h`, String)
print(script_output)
```  

## DownloadPDB.jl

```@repl
using MIToS
julia = Base.julia_cmd()
script_path = joinpath(pkgdir(MIToS), "scripts", "DownloadPDB.jl")
script_output = read(`$julia $script_path -h`, String)
print(script_output)
```  

## Distances.jl

```@repl
using MIToS
julia = Base.julia_cmd()
script_path = joinpath(pkgdir(MIToS), "scripts", "Distances.jl")
script_output = read(`$julia $script_path -h`, String)
print(script_output)
```  

## MSADescription.jl

```@repl
using MIToS
julia = Base.julia_cmd()
script_path = joinpath(pkgdir(MIToS), "scripts", "MSADescription.jl")
script_output = read(`$julia $script_path -h`, String)
print(script_output)
```  

## PercentIdentity.jl

```@repl
using MIToS
julia = Base.julia_cmd()
script_path = joinpath(pkgdir(MIToS), "scripts", "PercentIdentity.jl")
script_output = read(`$julia $script_path -h`, String)
print(script_output)
```  

## AlignedColumns.jl

```@repl
using MIToS
julia = Base.julia_cmd()
script_path = joinpath(pkgdir(MIToS), "scripts", "AlignedColumns.jl")
script_output = read(`$julia $script_path -h`, String)
print(script_output)
```  

## SplitStockholm.jl

```@repl
using MIToS
julia = Base.julia_cmd()
script_path = joinpath(pkgdir(MIToS), "scripts", "SplitStockholm.jl")
script_output = read(`$julia $script_path -h`, String)
print(script_output)
```

