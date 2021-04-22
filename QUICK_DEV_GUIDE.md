# Quick dev guide

If you are not very familiar with development in *Julia*, you can start with 
this simple approach.

1. Clone the repo from *GitHub* and enter the repo directory
2. Start *Julia REPL*
3. Change to *Pkg* mode in *Julia* (press `]`) and activate the environment for 
   the repo:
    ```julia
    pkg> activate .
    ```
4. Go back to normal REPL mode (press backspace) and load 
   [*Revise*](https://github.com/timholy/Revise.jl)
    ```julia
    julia> using Revise
    ```
5. Load *MIToS*
    ```julia
    julia> using MIToS
    ```
6. (optional) Check that *Revise* is tracking the correct files
    ```julia
    julia> Revise.watched_files
    ```

Edit the code, and the changes should be automatically loaded into the current 
session.


Happy coding!
