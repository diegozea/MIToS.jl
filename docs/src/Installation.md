
# Installation

First you need [**Julia 0.5.2**![](./assets/external-link.png)](https://julialang.org/downloads/oldreleases.html)
installed. A new version of MIToS with Julia 0.6 support is under development, we expect to
release it soon. The MIToS' stable version can be installed by typing on the Julia REPL:  

```julia
Pkg.add("MIToS")
```

If everything goes well with the installation, MIToS will be loaded without errors by typing:  

```julia
using MIToS
```

You can optionally do an exhaustive test of your installed version of MIToS with `Pkg.test` (it takes few seconds):  

```julia
Pkg.test("MIToS")
```

!!! note
    **Ways to run Julia**  

    *Option* | *Description*  
    ---:| ---  
    [Julia REPL![](./assets/external-link.png)](http://docs.julialang.org/en/stable/manual/getting-started/) | Built-in Julia command line. Start an Julia interactive session (REPL) by double-clicking the Julia executable or running `julia` from the system command line.
    [JuliaBox![](./assets/external-link.png)](https://www.juliabox.org/) | You can try Julia from your *web browser*. *No installation is required.*
    [IJulia![](./assets/external-link.png)](https://github.com/JuliaLang/IJulia.jl) | *Jupyter/IPython notebook* for Julia. It was used for generating the *this* documentation.
    [Juno![](./assets/external-link.png)](http://junolab.org/) | Integrated Development Environment (IDE).  



## Plots installation

Julia plotting capabilities are available through external packages. MIToS make use of
 *RecipesBase* to define plot recipes, which can be plotted using
 [Plots![](./assets/external-link.png)](https://juliaplots.github.io/) and different
 backends. You need to [install Plots![](./assets/external-link.png)](https://juliaplots.github.io/install/)
 to plot MIToS objects:  

```julia
Pkg.add("Plots")
```

And you also need to install at least one of the following backends:  

```julia
Pkg.add("GR") # Fast
Pkg.add("PlotlyJS") # Interactive
```

You need to load Plots in order to use the `plot` function. There is more information about
it in the [Plots documentation![](./assets/external-link.png)](https://juliaplots.github.io/).  

```julia
using Plots
```

To generate **graph** (network), **arc** and **chord** (circo) **plots**, you also need to
install and load [PlotRecipes![](./assets/external-link.png)](https://github.com/JuliaPlots/PlotRecipes.jl).  

```julia  
Pkg.add("PlotRecipes")

using PlotRecipes
```

## Scripts location

The MIToS’ scripts are located in the `MIToS/scripts` folder and can be runned from your
system command line. It’s possible to ask Julia for the location of the installed package
using `Pkg.dir`  


```julia
joinpath(Pkg.dir("MIToS"), "scripts")
```

You might want to add this folder into your `PATH` to easily access MIToS’ scripts.  

### How to add the script folder to `PATH` in Bash?

You can do it by adding the path of the MIToS script folder into the `~/.bashrc` file:

```julia
using MIToS
open(joinpath(homedir(), ".bashrc"), "r+") do fh
    path_to_scripts = joinpath(dirname(pathof(MIToS)), "..", "scripts")
    if all(line -> !occursin(path_to_scripts, line), eachline(fh))
        println(fh, "export PATH=\"\$PATH:", path_to_scripts, "\"")
    end
end
```
