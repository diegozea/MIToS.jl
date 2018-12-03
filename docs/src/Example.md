```@setup log
@info "Example"
```

# Example

In this simple demonstration, you will see how to calculate **ZBLMIp** (**Z** score of the
corrected **MIp** using BLOSUM62 pseudo frequencies) for a [Pfam![](./assets/external-link.png)](http://pfam.xfam.org/)
MSA from the [Julia REPL](@ref juliarepl) or using a
[MIToS script in the system command line](@ref commandline).  

## [MIToS in the Julia REPL](@id juliarepl)

If you load the `Pfam` module from `MIToS`, you will get access to a set of functions that
work with Pfam MSAs. In this case, we are going to use it for download a
[Stockholm![](./assets/external-link.png)](https://en.wikipedia.org/wiki/Stockholm_format)
MSA from the Pfam website and read it into Julia.  

```@setup juliarepl
using Plots
gr() # Just to avoid warnings in the output
```

```@example juliarepl
using MIToS.Pfam
pfam_file = downloadpfam("PF10660")
msa = read(pfam_file, Stockholm, generatemapping=true, useidcoordinates=true)
```

!!! note
    **Generation of sequence and column mappings**  
    The keyword argument `generatemapping` of `read` allows to generate sequence and column
    mappings for the MSA. *Column mapping* is the map between of each column on the MSA
    object and the column number in the file. *Sequence mappings* will use the start and
    end coordinates in the sequence ids for enumerate each residue in the sequence if
    `useidcoordinates` is `true`.  

You can plot this MSA and other MIToS’ objects using the [Plots![](./assets/external-link.png)](https://juliaplots.github.io/) package. The installation of *Plots* is described in the *Installation* section of this site:

```@example juliarepl
using Plots
gr()
plot(msa)
png("msa.png") # hide
nothing # hide
```  

![](msa.png)  

The `Information` module of `MIToS` has functions to calculate measures from the
[Information Theory![](./assets/external-link.png)](https://en.wikipedia.org/wiki/Information_theory),
such as Entropy and Mutual Information (MI), on a MSA. In this example, we will estimate
covariation between columns of the MSA with a corrected **MI** that use the BLOSUM62 matrix
for calculate pseudo frequencies (`BLMI`).  

```@example juliarepl
using MIToS.Information
ZBLMIp, BLMIp = BLMI(msa)
ZBLMIp # shows ZBLMIp scores
```

Once the *Plots* package is installed and loaded, you can use its capabilities to visualize
this results:

```@example juliarepl
heatmap(ZBLMIp, yflip=true, c=:grays)
png("blmi.png") # hide
nothing # hide
```  

![](blmi.png)  

```@setup juliarepl
rm(pfam_file) # clean up
```

## [MIToS in system command line](@id commandline)

Calculate ZBLMIp on the system shell is easy using the MIToS script called `BLMI.jl`. This
script reads a MSA file, and writes a file with the same base name of the input but with
the `.BLMI.csv` extension.  

```
BLMI.jl PF14972.stockholm.gz
```
