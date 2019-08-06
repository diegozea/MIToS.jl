# Quick dev guide

    If you are not very familiar with development in Julia, you can start with
    this simple approach.

1. Clone the repo from github and change to the repo directory
1. Start Julia REPL
1. Change to pkg mode in julia (press ']') and activate the enviroment for the repo.
    - pkg> activate .
1. Go back to normal REPL mode (press backspace) and load Revise
    julia> using Revise
1. Load MIToS
    julia> using MIToS
1. Check that Revise is tracking the correct files (optional)
    julia> Revise.watched_files

Edit the code and changes should be automatically loaded into the current session.


Happy coding.
