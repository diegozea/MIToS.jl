if [[ "$TRAVIS_OS_NAME" == "linux" ]]
then
    # COVERAGE
    julia -e 'cd(Pkg.dir("MIToS")); Pkg.add("Coverage"); using Coverage; Coveralls.submit(process_folder()); Codecov.submit(Codecov.process_folder())'
    # DOCS
    if [[ "$TRAVIS_BRANCH" != "master" ]]
    then
        julia -e 'Pkg.add("Documenter")'
        julia -e 'Pkg.add("Plots")' # Plots
        julia -e 'ENV["PYTHON"] = ""; Pkg.add("PyPlot"); Pkg.build("PyPlot")' # PyPlot
        julia -e 'Pkg.add("DataFrames"); Pkg.add("StatPlots")' # These are used in MSA.md
        julia -e 'cd(Pkg.dir("MIToS")); include(joinpath("docs", "make.jl"))'
    fi
fi
