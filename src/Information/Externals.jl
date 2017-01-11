# GaussDCA
# ========

"""
Wrapper function to `GaussDCA.gDCA`.
You need to install GaussDCA:
```julia
Pkg.clone("https://github.com/carlobaldassi/GaussDCA.jl")
```
Look into [GaussDCA.jl README](https://github.com/carlobaldassi/GaussDCA.jl) for further information.
If you use this wrapper, **please cite the GaussDCA publication and the package's doi**.

**GaussDCA Publication:**
Baldassi, Carlo, Marco Zamparo, Christoph Feinauer, Andrea Procaccini, Riccardo Zecchina, Martin Weigt, and Andrea Pagnani.
"Fast and accurate multivariate Gaussian modeling of protein families: predicting residue contacts and protein-interaction partners."
PloS one 9, no. 3 (2014): e92721.
"""
function gaussdca(msa; kargs...)
    plm = fill!(columnpairsmatrix(msa), NaN)
    msa_file = tempname() * ".fasta"
    write(msa_file, msa, FASTA)
    jdl_file = tempname() * ".jld"
    try
        bin_julia = ENV["_"]
        string_call = _create_string_call(msa_file, jdl_file; kargs...)
        run(`$bin_julia -e $string_call`)
        pairedvalues = JLD.load(jdl_file, "values")
        for (i,j,value) in pairedvalues
           plm[i,j] = value
        end
    finally
        isfile(msa_file) && rm(msa_file)
        isfile(jdl_file) && rm(jdl_file)
    end
    plm
end

function _create_string_call(msa_file::String, jdl_file::String; kargs...)
    if length(kargs) > 0
        for (k,v) in kargs
            @assert isa(v,Number) || isa(v,Symbol) "Argument values must be numbers or symbols"
        end
        str_kargs = "," * join([ isa(v,Symbol) ? "$k=:$v" : "$k=$v" for (k,v) in kargs],',')
    else
        str_kargs = ""
    end
    string(
    "using GaussDCA, JLD;",
    "values = GaussDCA.gDCA(\"$msa_file\"$str_kargs);",
    "JLD.save(\"$(jdl_file)\", \"values\", values)"
    )
end
