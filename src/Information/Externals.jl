# GaussDCA
# ========

"""
Wrapper function to `GaussDCA.gDCA`.
You need to install GaussDCA:
```julia
Pkg.clone("https://github.com/carlobaldassi/GaussDCA.jl")
```
And you need to load it before using this wrapper:
```julia
using GaussDCA
```
Look into [GaussDCA.jl README](https://github.com/carlobaldassi/GaussDCA.jl) for further information.
If you use this wrapper, **please cite the GaussDCA publication and the package's doi**.

**GaussDCA Publication:**
Baldassi, Carlo, Marco Zamparo, Christoph Feinauer, Andrea Procaccini, Riccardo Zecchina, Martin Weigt, and Andrea Pagnani.
"Fast and accurate multivariate Gaussian modeling of protein families: predicting residue contacts and protein-interaction partners."
PloS one 9, no. 3 (2014): e92721.
"""
function gaussdca(msa, args...; kargs...)
    filename = tempname()
    write(filename, msa, FASTA)
    try
        pairedvalues = Main.GaussDCA.gDCA(filename, args...; kargs...)
        plm = fill!(columnpairsmatrix(msa), NaN)
        for (i,j,value) in pairedvalues
           plm[i,j] = value
        end
        plm
    finally
        rm(filename)
    end
end

