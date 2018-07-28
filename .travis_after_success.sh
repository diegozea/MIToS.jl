if [[ "$TRAVIS_OS_NAME" == "linux" ]]
then
    # COVERAGE
    julia -e 'cd(Pkg.dir("MIToS")); Pkg.add("Coverage"); using Coverage; Coveralls.submit(process_folder()); Codecov.submit(Codecov.process_folder())'
    # DOCS
    # sudo python -m pip install --upgrade pip # PyPlot
    # sudo pip install --upgrade matplotlib # PyPlot
    # Work around matplotlib backend error leading to "matplotlib.pyplot could not be found by pyimport": can't just set env PYTHON="".
    julia -e 'ENV["PYTHON"]=""; Pkg.add("Conda"); using Conda; Conda.update(); Conda.add("matplotlib"); Pkg.add("PyCall"); Pkg.build("PyCall"); Pkg.add("PyPlot");'
    export LD_LIBRARY_PATH=$HOME/.julia/v0.6/Conda/deps/usr/lib
    julia -e 'using PyPlot' # PyPlot
    julia -e 'Pkg.add("Plots")' # Plots
    julia -e 'Pkg.add("PlotRecipes")' # PlotRecipes
    julia -e 'Pkg.checkout("PlotRecipes")' # PlotRecipes
    julia -e 'Pkg.add("DataFrames"); Pkg.add("StatPlots")' # These are used in MSA.md
    julia -e 'Pkg.add("Documenter")'
    echo "=============="
    echo "  DOCUMENTER  "
    echo "=============="
    julia -e 'cd(Pkg.dir("MIToS")); include(joinpath("docs", "make.jl"))'
fi
