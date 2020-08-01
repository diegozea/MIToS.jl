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
julia_path = joinpath(Base.Sys.BINDIR, "julia")
script_path = joinpath(dirname(pathof(MIToS)), "..", "scripts", "Buslje09.jl")
read(`$julia_path $script_path -h`, String) |> println
```  

## BLMI.jl

```@repl
using MIToS
julia_path = joinpath(Base.Sys.BINDIR, "julia")
script_path = joinpath(dirname(pathof(MIToS)), "..", "scripts", "BLMI.jl")
read(`$julia_path $script_path -h`, String) |> println
```  

## Conservation.jl

```@repl
using MIToS
julia_path = joinpath(Base.Sys.BINDIR, "julia")
script_path = joinpath(dirname(pathof(MIToS)), "..", "scripts", "Conservation.jl")
read(`$julia_path $script_path -h`, String) |> println
```  

## DownloadPDB.jl

```@repl
using MIToS
julia_path = joinpath(Base.Sys.BINDIR, "julia")
script_path = joinpath(dirname(pathof(MIToS)), "..", "scripts", "DownloadPDB.jl")
read(`$julia_path $script_path -h`, String) |> println
```  

## Distances.jl

```@repl
using MIToS
julia_path = joinpath(Base.Sys.BINDIR, "julia")
script_path = joinpath(dirname(pathof(MIToS)), "..", "scripts", "Distances.jl")
read(`$julia_path $script_path -h`, String) |> println
```  

## MSADescription.jl

```@repl
using MIToS
julia_path = joinpath(Base.Sys.BINDIR, "julia")
script_path = joinpath(dirname(pathof(MIToS)), "..", "scripts", "MSADescription.jl")
read(`$julia_path $script_path -h`, String) |> println
```  

## PercentIdentity.jl

```@repl
using MIToS
julia_path = joinpath(Base.Sys.BINDIR, "julia")
script_path = joinpath(dirname(pathof(MIToS)), "..", "scripts", "PercentIdentity.jl")
read(`$julia_path $script_path -h`, String) |> println
```  

## AlignedColumns.jl

```@repl
using MIToS
julia_path = joinpath(Base.Sys.BINDIR, "julia")
script_path = joinpath(dirname(pathof(MIToS)), "..", "scripts", "AlignedColumns.jl")
read(`$julia_path $script_path -h`, String) |> println
```  

## SplitStockholm.jl

```@repl
using MIToS
julia_path = joinpath(Base.Sys.BINDIR, "julia")
script_path = joinpath(dirname(pathof(MIToS)), "..", "scripts", "SplitStockholm.jl")
read(`$julia_path $script_path -h`, String) |> println
```

