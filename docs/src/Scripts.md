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
julia_path = joinpath(Base.JULIA_HOME, "julia")
script_path = joinpath(Pkg.dir("MIToS"), "scripts", "Buslje09.jl")
run(`$julia_path $script_path -h`)
```  

## BLMI.jl

```@repl
julia_path = joinpath(Base.JULIA_HOME, "julia")
script_path = joinpath(Pkg.dir("MIToS"), "scripts", "BLMI.jl")
run(`$julia_path $script_path -h`)
```  

## Conservation.jl

```@repl
julia_path = joinpath(Base.JULIA_HOME, "julia")
script_path = joinpath(Pkg.dir("MIToS"), "scripts", "Conservation.jl")
run(`$julia_path $script_path -h`)
```  

## DownloadPDB.jl

```@repl
julia_path = joinpath(Base.JULIA_HOME, "julia")
script_path = joinpath(Pkg.dir("MIToS"), "scripts", "DownloadPDB.jl")
run(`$julia_path $script_path -h`)
```  

## Distances.jl

```@repl
julia_path = joinpath(Base.JULIA_HOME, "julia")
script_path = joinpath(Pkg.dir("MIToS"), "scripts", "Distances.jl")
run(`$julia_path $script_path -h`)
```  

## MSADescription.jl

```@repl
julia_path = joinpath(Base.JULIA_HOME, "julia")
script_path = joinpath(Pkg.dir("MIToS"), "scripts", "MSADescription.jl")
run(`$julia_path $script_path -h`)
```  

## PercentIdentity.jl

```@repl
julia_path = joinpath(Base.JULIA_HOME, "julia")
script_path = joinpath(Pkg.dir("MIToS"), "scripts", "PercentIdentity.jl")
run(`$julia_path $script_path -h`)
```  

## AlignedColumns.jl

```@repl
julia_path = joinpath(Base.JULIA_HOME, "julia")
script_path = joinpath(Pkg.dir("MIToS"), "scripts", "AlignedColumns.jl")
run(`$julia_path $script_path -h`)
```  

## SplitStockholm.jl

```@repl
julia_path = joinpath(Base.JULIA_HOME, "julia")
script_path = joinpath(Pkg.dir("MIToS"), "scripts", "SplitStockholm.jl")
run(`$julia_path $script_path -h`)
```  
