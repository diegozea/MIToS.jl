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

To update MIToS to the latest version, you can run:

```julia
using Pkg
Pkg.update("MIToS")
```

!!! tip "Ways to run Julia"
    

  - **[Julia REPL ![](./assets/external-link.png)](https://docs.julialang.org/en/v1/stdlib/REPL/):** Built-in Julia command line. Start a Julia interactive session (REPL) by double-clicking the Julia executable or running `julia` from the system command line.
  - **[IJulia ![](./assets/external-link.png)](https://github.com/JuliaLang/IJulia.jl):** *Jupyter/IPython notebook* for Julia.
  - **[Pluto ![](./assets/external-link.png)](https://github.com/fonsp/Pluto.jl):** A simple reactive notebook for Julia.
  - **[VS Code Extension for Julia ![](./assets/external-link.png)](https://www.julia-vscode.org/):** The Julia's Integrated Development Environment (IDE).

!!! info "Running the test suite"
    
    **Optionally**, you can run the test suite to ensure everything works as expected.
    The test suite is extensive and can take several minutes to run. It is the same test
    suite used for MIToS' continuous integration (CI), so everything should pass. To run
    the test suite, execute `using Pkg; Pkg.test("MIToS")` in the Julia REPL.

## Plots installation

Julia plotting capabilities are available through external packages. MIToS makes use of
*RecipesBase* to define plot recipes, which can be plotted using
[Plots![](./assets/external-link.png)](http://docs.juliaplots.org/latest/) and its different
backends. You need to [install Plots![](./assets/external-link.png)](http://docs.juliaplots.org/latest/install/)
to plot MIToS objects:

```julia
using Pkg
Pkg.add("Plots")
```

Once it is installed, you need to load Plots in order to use the `plot` function. There is
more information about it in the [Plots documentation![](./assets/external-link.png)](http://docs.juliaplots.org/latest/).

```julia
using Plots
```

To generate **graph** (network), **arc** and **chord** (circo) **plots**, you also need to
install and load [GraphRecipes![](./assets/external-link.png)](https://github.com/JuliaPlots/GraphRecipes.jl).

```julia
using Pkg
Pkg.add("GraphRecipes")

using GraphRecipes
```

You can look for examples in the [GraphRecipes documentation![](./assets/external-link.png)](https://docs.juliaplots.org/stable/GraphRecipes/examples/).
