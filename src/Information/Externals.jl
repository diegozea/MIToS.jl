# GaussDCA
# ========

"""
Wrapper function to `GaussDCA.gDCA`.
You need to install GaussDCA:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/carlobaldassi/GaussDCA.jl"))
```
Look into [GaussDCA.jl README](https://github.com/carlobaldassi/GaussDCA.jl) for further information.
If you use this wrapper, **please cite the GaussDCA publication and the package's doi**.

It's possible to indicate the path to the julia binary where GaussDCA is installed.
However, it's recommended to use the same version where MIToS is installed. That is
because this function use `serialize`/`deserialize` to transfer data between the processes.

**GaussDCA Publication:**
Baldassi, Carlo, Marco Zamparo, Christoph Feinauer, Andrea Procaccini, Riccardo Zecchina, Martin Weigt, and Andrea Pagnani.
"Fast and accurate multivariate Gaussian modeling of protein families: predicting residue contacts and protein-interaction partners."
PloS one 9, no. 3 (2014): e92721.
"""
function gaussdca(msa; juliapath::String=joinpath(Sys.BINDIR,Base.julia_exename()), kargs...)
    base_name = tempname()
    if Sys.iswindows()
        base_name = escape_string(base_name)
        juliapath = escape_string(juliapath)
        if !endswith(juliapath, ".exe")
            juliapath = juliapath * ".exe"
        end
    end
    script_file = base_name * ".jl"
    msa_file = base_name * ".fasta"
    jdl_file = base_name * ".jls"
    plm = fill!(columnpairsmatrix(msa), NaN)
    write(msa_file, msa, FASTA)
    try
        _create_script(script_file, msa_file, jdl_file; kargs...)
        run(`$juliapath $script_file`)
        pairedvalues = open(deserialize, jdl_file, "r")
        for (i,j,value) in pairedvalues
           plm[i,j] = value
        end
    finally
        isfile(script_file) && rm(script_file)
        isfile(msa_file) && rm(msa_file)
        isfile(jdl_file) && rm(jdl_file)
    end
    plm
end

function _create_script(script_file::String, msa_file::String, jdl_file::String; kargs...)
    if length(kargs) > 0
        for (k,v) in kargs
            @assert isa(v,Number) || isa(v,Symbol) "Argument values must be numbers or symbols"
        end
        str_kargs = "," * join([ isa(v,Symbol) ? "$k=:$v" : "$k=$v" for (k,v) in kargs],',')
    else
        str_kargs = ""
    end
    open(script_file, "w") do fh
        println(fh, "using Serialization;")
        println(fh, "using GaussDCA;")
        println(fh, "values = GaussDCA.gDCA(\"$msa_file\"$str_kargs);")
        println(fh, "open(fh -> serialize(fh, values), \"$(jdl_file)\", \"w\")")
    end
end
