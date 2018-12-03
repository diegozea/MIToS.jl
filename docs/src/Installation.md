```@setup log
@info "Installation docs"
```

# Installation

First you need to install [**Julia.**![](./assets/external-link.png)](https://julialang.org/downloads/)
MIToS' stable version can be installed by typing on the Julia REPL:  

```julia
using Pkg
Pkg.add("MIToS")
```

If everything goes well with the installation, MIToS will be loaded without errors by typing:  

```julia
using MIToS
```

You can optionally do an exhaustive test of your installed version of MIToS with `Pkg.test` (it takes few minutes):  

```julia
using Pkg
Pkg.test("MIToS")
```

!!! note
    **Ways to run Julia**  

    *Option* | *Description*  
    ---:| ---  
    [Julia REPL![](./assets/external-link.png)](https://docs.julialang.org/en/v1/stdlib/REPL/) | Built-in Julia command line. Start an Julia interactive session (REPL) by double-clicking the Julia executable or running `julia` from the system command line.
    [JuliaBox![](./assets/external-link.png)](https://juliabox.com/) | You can try Julia from your *web browser*. *No installation is required.*
    [IJulia![](./assets/external-link.png)](https://github.com/JuliaLang/IJulia.jl) | *Jupyter/IPython notebook* for Julia.
    [Juno![](./assets/external-link.png)](http://junolab.org/) | Integrated Development Environment (IDE).  



## Plots installation

Julia plotting capabilities are available through external packages. MIToS makes use of
 *RecipesBase* to define plot recipes, which can be plotted using
 [Plots![](./assets/external-link.png)](http://docs.juliaplots.org/latest/) and different
 backends. You need to [install Plots![](./assets/external-link.png)](http://docs.juliaplots.org/latest/install/)
 to plot MIToS objects:  

```julia
using Pkg
Pkg.add("Plots")
```

And you also need to install at least one of the following backends:  

```julia
using Pkg
Pkg.add("GR") # Fast
Pkg.add("PlotlyJS") # Interactive
```

You need to load Plots in order to use the `plot` function. There is more information about
it in the [Plots documentation![](./assets/external-link.png)](http://docs.juliaplots.org/latest/).  

```julia
using Plots
```

To generate **graph** (network), **arc** and **chord** (circo) **plots**, you also need to
install and load [GraphRecipes![](./assets/external-link.png)](https://github.com/JuliaPlots/GraphRecipes.jl).  

```julia  
Pkg.add("GraphRecipes")

using GraphRecipes
```

## Scripts location

The MIToS’ scripts are located in the `MIToS/scripts` folder and can be runned from your
system command line. It’s possible to ask Julia for the location of the installed package
using:


```julia
import MIToS
joinpath(splitdir(dirname(pathof(MIToS)))[1], "scripts")
```

You might want to add this folder into your `PATH` to easily access MIToS’ scripts.  
For example, in **bash** you can do it by adding the path of the MIToS script folder
into the `~/.bashrc` file. The `println` output shows the line to add to that file:

```julia
import MIToS
println("export PATH=\"\$PATH:", joinpath(splitdir(dirname(pathof(MIToS)))[1], "scripts"), "\"")
```
