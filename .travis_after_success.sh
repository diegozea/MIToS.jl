if [[ "$TRAVIS_OS_NAME" == "linux" ]]
then
    # COVERAGE
    julia -e 'cd(Pkg.dir("MIToS")); Pkg.add("Coverage"); using Coverage; Coveralls.submit(process_folder()); Codecov.submit(Codecov.process_folder())'
    # DOCS
    sudo python -m pip install --upgrade pip # PyPlot
    sudo pip install --upgrade matplotlib # PyPlot
    julia -e 'Pkg.add("PyPlot"); Pkg.build("PyPlot")' # PyPlot
    julia -e 'Pkg.add("Plots")' # Plots
    julia -e 'Pkg.add("PlotRecipes")' # PlotRecipes
    julia -e 'Pkg.checkout("PlotRecipes")' # PlotRecipes
    julia -e 'Pkg.add("DataFrames"); Pkg.add("StatPlots")' # These are used in MSA.md
    julia -e 'Pkg.add("Documenter")'
    echo "DOCS"
    julia -e 'cd(Pkg.dir("MIToS")); include(joinpath("docs", "make.jl"))'
fi
