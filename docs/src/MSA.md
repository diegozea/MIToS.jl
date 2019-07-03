```@setup log
@info "MSA docs"
```

# [MSA](@id Module-MSA)

The MSA module of MIToS has utilities for working with Multiple Sequence Alignments of
protein Sequences (MSA).  

```julia
using MIToS.MSA # to load the MSA module
```

## Features

- [**Read**](@ref Reading-MSA-files) and [**write**](@ref Writing-MSA-files) MSAs in `Stockholm`, `FASTA`, `PIR` or `Raw` format.
- Handle [**MSA annotations**](@ref MSA-Annotations).
- [**Edit the MSA**](@ref Editing-your-MSA), e.g. delete columns or sequences, change sequence order, shuffling...
- [**Keep track of positions**](@ref Column-and-sequence-mappings) and annotations after modifications on the MSA.
- [**Describe an MSA**](@ref Describing-your-MSA), e.g. mean percent identity, sequence coverage, gap percentage...
- [**Sequence clustering**](@ref Sequence-clustering) with Hobohm I.

## Contents

```@contents
Pages = ["MSA.md"]
Depth = 4
```  

## [MSA IO](@id MSA-IO)

### [Reading MSA files](@id Reading-MSA-files)

The main function for reading MSA files in MIToS is `read` and it is defined in the `Utils`
module. This function takes a filename/path as a first argument followed by other
arguments. It opens the file and uses the arguments to call the `parse` function.
`read` decides how to open the file, using the prefixes (e.g. https) and suffixes
(i.e. extensions) of the file name, while `parse` does the actual parsing of
the file. You can `read` **gzipped files** if they have the `.gz` extension and
also urls pointing to a **web file**.  
The second argument of `read` and `parse` is the file `FileFormat`. The supported MSA formats
at the moment are `Stockholm`, `FASTA`, `PIR` (NBRF) and `Raw`.  
For example, reading with MIToS the full Stockholm MSA of the family PF07388 using the Pfam
RESTful interface will be:  

```@example msa_read
using MIToS.MSA

read("http://pfam.xfam.org/family/PF07388/alignment/full", Stockholm)
```

The third (and optional) argument of `read` and `parse` is the output MSA type:  

- `Matrix{Residue}` : It only contains the aligned sequences.  
- `MultipleSequenceAlignment` : It contains the aligned sequences and their
names/identifiers.  
- `AnnotatedMultipleSequenceAlignment` : It's the richest MIToS' MSA format and it's the
default. It includes the aligned sequences, their names and the MSA annotations.  

Example of `Matrix{Residue}` output using a `Stockholm` file as input:

```@example msa_read
read("http://pfam.xfam.org/family/PF07388/alignment/full", Stockholm, Matrix{Residue})
```

Because `read` calls `parse`, you should look into the documentation of `parse`
to know the available keyword arguments. The optional keyword arguments of
those functions are:  

- `generatemapping` : If `generatemapping` is `true` (default: `false`), sequences and
columns mappings are generated and saved in the MSA annotations. **The default is `false` to
not overwrite mappings by mistake when you read an annotated MSA file saved with MIToS.**  
- `useidcoordinates` : If `useidcoordinates` is `true` (default: `false`) and the names
have the form *seqname/start-end*, MIToS uses this coordinates to generate sequence
mappings. This is safe and useful with unmodified Pfam MSAs. **Do not use it when reading an
MSA saved with MIToS. MIToS deletes unaligned insert columns, therefore disrupts sequences
that have them.**  
- `deletefullgaps` : Given that lowercase characters and dots are converted to gaps,
unaligned insert columns in the MSA (derived from a HMM profile) are converted into full
gap columns. `deletefullgaps` is `true` by default, deleting full gaps columns and
therefore insert columns.  

!!! note  
    **If you want to keep the insert columns...**  Use the keyword argument `keepinserts`
    to `true` in `read`/`parse`. This only works with an `AnnotatedMultipleSequenceAlignment`
    output. A column annotation (`"Aligned"`) is stored in the annotations, where insert
    columns are marked with `0` and aligned columns with `1`.  

When `read` returns an `AnnotatedMultipleSequenceAlignment`, it uses the MSA `Annotations`
to keep track of performed modifications. To access these notes, use `printmodifications`:  

```@example msa_read
msa = read("http://pfam.xfam.org/family/PF01565/alignment/full", Stockholm)

printmodifications(msa)
```  

### [Writing MSA files](@id Writing-MSA-files)

Julia REPL shows MSAs as Matrices. If you want to print them in another format, you should
use the `print` function with an MSA object as first argument and the `FileFormat` `FASTA`,
`Stockholm`, `PIR` or `Raw` as second argument.  

```@example msa_write
using MIToS.MSA

msa = read("http://pfam.xfam.org/family/PF16996/alignment/full", Stockholm) # reads a Stockholm MSA file

print(msa, FASTA) # prints msa in FASTA format
```  

To save an MSA object to a file, use the `write` function. This function takes a filename
as a first argument. If the filename ends with `.gz`, the output will be a compressed
(gzipped) file. The next two arguments of `write` are passed to `print`, so `write` behaves
as `print`.  

```@example msa_write
write("msa.gz", msa, FASTA) # writes msa in FASTA format in a gzipped file
```  

## [MSA Annotations](@id MSA-Annotations)

MSA annotations are based on the Stockholm format mark-ups. There are four types of
annotations stored as dictionaries. All the annotations have a feature name as part of the
key, which should be a single "word" (without spaces) and less than 50 characters long.  

- **File annotations** : The annotations can contain either file or MSA information. They
have feature names as keys and the values are strings (free text). Lines starting with
`#=GF` in Stockholm format.  
- **Column annotations** : They have feature names as keys and strings with exactly 1 char
per column as values. Lines starting with `#=GC` in Stockholm format.  
- **Sequence annotations** : The keys are tuples with the sequence name and the feature
name. The values are free text (strings). Lines starting with `#=GS` in Stockholm format.
Annotations in the `PIR`/NBRF format are also stored as sequence annotations. In
particular, we use the names `"Type"` and `"Title"` to name the sequence type in the
identifier line and the first comment line before the sequence in PIR files, respectively.  
- **Residue annotations** : The keys are tuples with the sequence name and the feature
name. The values are strings with exactly 1 char per column/residues. `#=GR` lines in
Stockholm format.  

Julia REPL shows the `Annotations` type as they are represented in the [Stockholm format![](./assets/external-link.png)](https://en.wikipedia.org/wiki/Stockholm_format).
You can get the `Annotations` inside an annotated MSA or sequence using the `annotations`
function.  

```@example msa_write
annotations(msa)
```  

Particular annotations can be accessed using the functions `getannot...`. These functions
take the MSA/sequence as first argument and the feature name of the desired annotation as
the last. In the case of `getannotsequence` and `getannotresidue`, the second argument
should be the sequence name.  

```@example msa_write
getannotsequence(msa, "A0A139NPI6_9STRE/5-59", "AC") # ("A0A139NPI6_9STRE/5-59", "AC") is the key in the dictionary
```  

If you want to add new annotations, you should use the `setannot…!` functions. These
functions have the same arguments that `getannot...` functions except for an
extra argument used to indicate the new annotation value.  

```@example msa_write
setannotsequence!(msa, "A0A139NPI6_9STRE/5-59", "New_Feature_Name", "New_Annotation")
```  

A `getannot...` function without the key (last arguments), returns the particular
annotation dictionary. As you can see, the new sequence annotation is now part of our
MSA annotations.  

```@example msa_write
getannotsequence(msa)
```

## [Editing your MSA](@id Editing-your-MSA)

MIToS offers functions to edit your MSA. Because these functions modify the msa, their
names end with a bang `!`, following the Julia convention. Some of these functions have an
`annotate` keyword argument (in general, it's `true` by default) to indicate if the
modification should be recorded in the MSA/sequence annotations.  

One common task is to delete sequences or columns of the MSA. This could be done using the
functions `filtersequences!` and `filtercolumns!`. These functions take the MSA or sequence
(if it's possible) as first argument and a `BitVector` or `Vector{Bool}` mask as second
argument. It deletes all the sequences or columns where the mask is `false`. These functions
are also defined for `Annotations`, this allows to automatically update (modify) the
annotations (and therefore, sequence and column mappings) in the MSA.  

This two deleting operations are used in the second and third mutating
functions of the following list:  

- `setreference!` : Sets one of the sequences as the first sequence of the MSA (query or
reference sequence).  
- `adjustreference!` : Deletes columns with gaps in the first sequence of the MSA
(reference).  
- `gapstrip!` : This function first calls `adjustreference!`, then deletes sequences with
low (user defined) MSA coverage and finally, columns with user defined % of gaps.  

Also, there are several available funtions `shuffle_…!`. These functions are useful to
generate random alignments. The `Information` module of `MIToS` uses them to calculate the
Z scores of MI values.  

#### [Example: Deleting sequences](@id Example:-Deleting-sequences)

For example, if you want to keep only the proteins from *Actinobacteria* you can delete
all the sequences that don't have `_9ACTN` in their UniProt entry names:  

```@example msa_edit
using MIToS.MSA

msa = read("http://pfam.xfam.org/family/PF07388/alignment/full", Stockholm)

sequencenames(msa) # the function sequencenames returns the sequence names in the MSA
```  

```@example msa_edit
mask = map(x -> occursin(r"_9ACTN", x), sequencenames(msa)) # an element of mask is true if "_9ACTN" is in the name
```  

```@example msa_edit
filtersequences!(msa, mask) # deletes all the sequences where mask is false

sequencenames(msa)
```

#### [Example: Exporting a MSA for freecontact (part I)](@id Example:-Exporting-a-MSA-for-freecontact-(part-I))

The most simple input for the command line tool [freecontact![](./assets/external-link.png)](https://rostlab.org/owiki/index.php/FreeContact)
(if you don't want to set `--mincontsep`) is a `Raw` MSA file with a reference sequence
without insertions or gaps. This is easy to get with MIToS using `read` (deletes the insert
columns), `setreference!` (to choose a reference), `adjustreference!` (to delete columns
with gaps in the reference) and `write` (to save it in `Raw` format) functions.  

```@repl
using MIToS.MSA
msa = read("http://pfam.xfam.org/family/PF02476/alignment/full", Stockholm)
msa_coverage = coverage(msa)
maxcoverage, maxindex = findmax(msa_coverage) # chooses the sequence with more coverage of the MSA
setreference!(msa, maxindex[1])
adjustreference!(msa)
write("tofreecontact.msa", msa, Raw)
print(read("tofreecontact.msa", String)) # It displays the contents of the output file
```

## [Column and sequence mappings](@id Column-and-sequence-mappings)

Inserts in a Stockholm MSA allow to access the full fragment of the aligned sequences.
Using this, combined with the sequence names that contain coordinates used in Pfam, you
can know what is the UniProt residue number of each residue in the MSA.   

```julia
"PROT_SPECI/3-15 .....insertALIGNED"
#                     3456789111111
#                            012345
```

MIToS `read` and `parse` functions delete the insert columns, but they do the mapping
between each residue and its residue number before deleting insert columns when `generatemapping` is
`true`. If you don't set `useidcoordinates` to `true`, the residue first `i` residue will
be 1 instead of 3 in the previous example.  

```@example msa_mapping
using MIToS.MSA

msa = parse("PROT_SPECI/3-15 .....insertALIGNED", Stockholm, generatemapping=true, useidcoordinates=true)
```  

MIToS also keeps the column number of the input MSA and its total number of columns. All
this data is stored in the MSA annotations using the `SeqMap`, `ColMap` and `NCol` feature
names.  

```@example msa_mapping
annotations(msa)
```

To have an easy access to mapping data, MIToS provides the `getsequencemapping` and
`getcolumnmapping` functions.  

```@example msa_mapping
getsequencemapping(msa, "PROT_SPECI/3-15")
```  

```@example msa_mapping
getcolumnmapping(msa)
```  

#### [Example: Exporting a MSA for freecontact (part II)](@id Example:-Exporting-a-MSA-for-freecontact-(part-II))

If we want to use the `--mincontsep` argument of `freecontact` to calculate scores between
distant residues, we will need to add a header to the MSA. This header should contains the
residue number of the first residue of the sequence and the full fragment of that sequence
(with the inserts). This data is used by FreeContact to calculate the residue number of
each residue in the reference sequence.  
We are going to use MIToS mapping data to create this header, so we read the MSA with
`generatemapping` and `useidcoordinates` set to `true`.  

```@example freecontact_ii
using MIToS.MSA

msa = read( "http://pfam.xfam.org/family/PF02476/alignment/full", Stockholm,
            generatemapping=true, useidcoordinates=true)
```  

Here, we are going to choose the sequence with more coverage of the MSA as our reference
sequence.  

```@example freecontact_ii
msa_coverage = coverage(msa)
maxcoverage, maxindex = findmax(msa_coverage)
setreference!(msa, maxindex[1])
adjustreference!(msa)
```

MIToS deletes the residues in insert columns, so we are going to use the
sequence mapping to generate the whole fragment of the reference sequence
(filling the missing regions with `'x'`).  

```@example freecontact_ii
seqmap = getsequencemapping(msa, 1) # seqmap will be a vector with the residue numbers of the first sequence (reference)

seq = collect( stringsequence(msa, 1) ) # seq will be a Vector of Chars with the reference sequence

sequence = map(seqmap[1]:seqmap[end]) do seqpos # for each position in the whole fragment
    if seqpos in seqmap                         # if that position is in the MSA
        popfirst!(seq)                          # the residue is taken from seq
    else                                        # otherwise
        'x'                                     # 'x' is included
    end
end

sequence = join(sequence) # join the Chars on the Vector to create a string
```

Once we have the whole fragment of the sequence, we create the file and write the header in
the required format (as in the man page of freecontact).  

```@example freecontact_ii
open("tofreecontact.msa", "w") do fh

    println(fh, "# querystart=", seqmap[1])

    println(fh, "# query=", sequence )

end
```

As last (optional) argument, `write` takes the mode in which is opened the file. We use
`"a"` here to append the MSA to the header.  

```@example freecontact_ii
write("tofreecontact.msa", msa, Raw, "a")
```  

```@example freecontact_ii
print(read("tofreecontact.msa", String)) # It displays the contents of the output file
```

## [Get sequences from a MSA](@id Get-sequences-from-a-MSA)

It's possible to index the MSA as any other matrix to get an aligned sequence. This will be
return a `Array` of `Residue`s without annotations but keeping names/identifiers.  

```@example msa_indexing
using MIToS.MSA

msa = read( "http://pfam.xfam.org/family/PF16996/alignment/full", Stockholm,
            generatemapping=true, useidcoordinates=true)
```  

```@example msa_indexing
msa[2,:] # second sequence of the MSA, it keeps column names
```  

```@example msa_indexing
msa[2:2,:] # Using the range 2:2 to select the second sequence, keeping also the sequence name
```

If you want to obtain the aligned sequence with its name and annotations (and therefore
sequence and column mappings), you should use the function `getsequence`. This function
returns an `AlignedSequence` with the sequence name from a `MultipleSequenceAlignment` or
an `AnnotatedAlignedSequence`, that also contains annotations, from an
`AnnotatedMultipleSequenceAlignment`.  

```@example msa_indexing
secondsequence = getsequence(msa, 2)
```  

```@example msa_indexing
annotations(secondsequence)
```  

Use `stringsequence` if you want to get the sequence as a string.  

```@example msa_indexing
stringsequence(msa, 2)
```

Because matrices are stored columnwise in Julia, you will find useful the
`getresiduesequences` function when you need to heavily operate over sequences.  

```@example msa_indexing
getresiduesequences(msa)
```

## [Describing your MSA](@id Describing-your-MSA)

The MSA module has a number of functions to gain insight about your MSA. Using `MIToS.MSA`,
one can easily ask for...  

- The **number of columns and sequences** with the `ncolumns` and `nsequences` functions.  
- The fraction of columns with residues (**coverage**) for each sequence making use of the
`coverage` method.  
- The **fraction or percentage of gaps/residues** using with the functions `gapfraction`,
`residuefraction` and `columngapfraction`.  
- The **percentage of identity** (PID) between each sequence of the MSA or its mean value
with `percentidentity` and `meanpercentidentity`.  

The percentage identity between two aligned sequences is a common measure of sequence
similarity and is used by the `hobohmI` method to estimate and reduce MSA redundancy.
**MIToS functions to calculate percent identity don't align the sequences, they need
already aligned sequences.** Full gaps columns don't count to the alignment length.  

```@example msa_describe
using MIToS.MSA

msa = permutedims(
        hcat(   res"--GGG-",      # res"..." uses the @res_str macro to create a (column) Vector{Residue}
                res"---GGG" ), (2,1))
#        identities 000110 sum 2
#  aligned residues 001111 sum 4
```  

```@example msa_describe
percentidentity(msa[1,:], msa[2,:]) # 2 / 4
```  

To quickly calculate if the percentage of identity is greater than a determined value, use
that threshold as third argument. `percentidentity(seqa, seqb, pid)` is a lot more faster
than `percentidentity(seqa, seqb) >= pid`.  

```@example msa_describe
percentidentity(msa[1,:], msa[2,:], 62) # 50% >= 62%
```

#### [Example: Plotting gap percentage per column and coverage per sequence](@id Example:-Plotting-gap-percentage-per-column-and-coverage-per-sequence)

The `gapfraction` and `coverage` functions return a vector of numbers between `0.0` and
`1.0` (fraction of...). Sometime it's useful to plot this data to quickly understand the
MSA structure. In this example, we are going to use the [Plots![](./assets/external-link.png)](http://plots.readthedocs.org/en/latest/)
package for plotting, with the [GR![](./assets/external-link.png)](https://github.com/jheinen/GR.jl)
backend, but you are free to use any of the Julia plotting libraries.  

```@setup msa_plots
@info "MSA: Plots"
using Plots
gr() # Hide possible warnings
```

```@example msa_plots
using MIToS.MSA

msa = read("http://pfam.xfam.org/family/PF09776/alignment/full", Stockholm)

using Plots

gr(size=(600,300))

plot(   1:ncolumns(msa), # x is a range from 1 to the number of columns
        vec(columngapfraction(msa)) .* 100.0, # y is a Vector{Float64} with the percentage of gaps of each column
        linetype = :line,
        ylabel = "gaps [%]",
        xlabel = "columns",
        legend=false)

png("msa_gaps.png") # hide
nothing # hide
```  

![](msa_gaps.png)  

```@example msa_plots
plot(   1:nsequences(msa), # x is a range from 1 to the number of sequences
        vec(coverage(msa)) .* 100, # y is a Vector{Float64} with the coverage of each sequence
        linetype = :line,
        ylabel = "coverage [%]",
        xlabel = "sequences",
        legend=false)

png("msa_coverage.png") # hide
nothing # hide
```  

![](msa_coverage.png)  

```@example msa_plots
plot(msa)
png("msa_msa.png") # hide
nothing # hide
```  

![](msa_msa.png)  

#### [Example: Filter sequences per coverage and columns per gap fraction](@id Example:-Filter-sequences-per-coverage-and-columns-per-gap-fraction)

Taking advantage of the `filter...!` functions and the `coverage` and `columngapfraction`
functions, it's possible to delete short sequences or columns with a lot of gaps.  

```@example msa_plots
println("\tsequences\tcolumns")
println( "Before:\t", nsequences(msa), "\t\t", ncolumns(msa)  )
# delete sequences with less than 90% coverage of the MSA length:
filtersequences!(msa, coverage(msa) .>= 0.9)
# delete columns with more than 10% of gaps:
filtercolumns!(msa, columngapfraction(msa) .<= 0.1)
println( "After:\t", nsequences(msa), "\t\t",  ncolumns(msa)  )
```

```@example msa_plots
histogram(  vec(columngapfraction(msa)),
            # Using vec() to get a Vector{Float64} with the fraction of gaps of each column
            xlabel = "gap fraction in [0,1]", bins = 20, legend = false)
png("msa_hist_gaps.png") # hide
nothing # hide
```  

![](msa_hist_gaps.png)  

```@example msa_plots
histogram(  vec(coverage(msa) .* 100.0), #  Column with the coverage of each sequence
            xlabel = "coverage [%]", legend=false)
png("msa_hist_coverage.png") # hide
nothing # hide
```  

![](msa_hist_coverage.png)  


#### [Example: Plotting the percentage of identity between sequences](@id Example:-Plotting-the-percentage-of-identity-between-sequences)

The distribution of the percentage of identity between every pair of sequences in an MSA,
gives an idea of the MSA diversity. In this example, we are  using `percentidentity`
over an MSA to get those identity values.  

```@example msa_pid
using MIToS.MSA
msa = read("http://pfam.xfam.org/family/PF09776/alignment/full", Stockholm)
pid = percentidentity(msa)
nothing # hide
```

MIToS stores the matrix of percentage of identity between the aligned sequences as a
PairwiseListMatrix from the [PairwiseListMatrices![](./assets/external-link.png)](http://diegozea.github.io/PairwiseListMatrices.jl/)
package. This matrix type saves RAM, allowing the storage of  big matrices. In this
example, we use the `to_table` function of *PairwiseListMatrices* to convert the matrix
into a table with indices.  

```@example msa_pid
using PairwiseListMatrices

pidtable = to_table(pid, diagonal=false)
```

The function `quantile` gives a quick idea of the percentage identity distribution of the MSA.  

```@example msa_pid
using Statistics

quantile(convert(Vector{Float64}, pidtable[:,3]), [0.00, 0.25, 0.50, 0.75, 1.00])
```

The function `meanpercentidentity` gives the mean value of the percent identity
distribution for MSA with less than 300 sequences, or a quick estimate (mean PID in a
random sample of sequence pairs) otherwise unless you set `exact` to `true`.  

```@example msa_pid
meanpercentidentity(msa)
```

One can easily plot that matrix and its distribution using the `heatmap` and `histogram`
functions of the [Plots![](./assets/external-link.png)](https://github.com/tbreloff/Plots.jl)
package.  

```@setup msa_pid
@info "MSA: PID"
using Plots
gr() # Hide possible warnings
```

```@example msa_pid
using Plots
gr()
heatmap(convert(Matrix, pid), yflip=true, ratio=:equal)
png("msa_heatmap_pid.png") # hide
nothing # hide
```  

![](msa_heatmap_pid.png)  

```@example msa_pid
histogram(pidtable[:,3], xlabel ="Percentage of identity", legend=false)
png("msa_hist_pid.png") # hide
nothing # hide
```  

![](msa_hist_pid.png)  

## [Sequence clustering](@id Sequence-clustering)  

The `MSA` module allows to clusterize sequences in an MSA. The `hobohmI` function takes as
input an MSA followed by an identity threshold value, and returns a `Clusters` type
with the result of a [Hobohm I![](./assets/external-link.png)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2142204/)
sequence clustering. The Hobohm I algorithm will add a sequence to an existing cluster, if
the percentage of identity is equal or greater than the threshold.  
The `Clusters` is sub-type of `ClusteringResult` from the [Clustering.jl![](./assets/external-link.png)](http://clusteringjl.readthedocs.org/en/latest/index.html)
package. One advantage of use a sub-type of `ClusteringResult`is that you are able to use
any method defined on `Clustering.jl` like `varinfo` (Variation of Information) for example.
Also, you can use any clustering algorithm included in *Clustering.jl*, and convert its
result to an `Clusters` object to use it with MIToS.  
`MSA` defines the functions `nclusters` to get the resulting number of clusters, `counts`
to get the number of sequences on each cluster and `assignments` to get the cluster number
of each sequence. The most important method is `getweight`, which returns the weight of
each sequence. This method is used in the `Information` module of MIToS to reduce redundancy.  

#### [Example: Reducing redundancy of a MSA](@id Example:-Reducing-redundancy-of-a-MSA)

MSAs can suffer from an unnatural sequence redundancy and a high number of protein
fragments. In this example, we are using a sequence clustering to make a non-redundant set
of representative sequences. We are going to use the function `hobohmI` to perform the clustering with the
Hobohm I algorithm at 62% identity.  

```@setup msa_clusters
@info "MSA: Clusters"
using Plots
using StatsPlots
using DataFrames
gr() # Hide possible warnings
```

```@example msa_clusters
using MIToS.MSA

msa = read("http://pfam.xfam.org/family/PF09776/alignment/full", Stockholm)

println("This MSA has ", nsequences(msa), " sequences...")
```

```@example msa_clusters
clusters = hobohmI(msa, 62)
```

```@example msa_clusters
println("...but has only ", nclusters(clusters), " sequence clusters after a clustering at 62% identity.")
```

```@example msa_clusters
using Plots
gr()

plot(msa)
png("msa_clusters_i.png") # hide
nothing # hide
```  

![](msa_clusters_i.png)  

We are going to use the [DataFrames![](./assets/external-link.png)](http://dataframesjl.readthedocs.org/en/latest/)
package to easily select the sequence with the highest coverage of each cluster.  

```@example msa_clusters
using DataFrames

df = DataFrame( seqnum = 1:nsequences(msa),
                seqname = sequencenames(msa),
                cluster = assignments(clusters), # the cluster number/index of each sequence
                coverage = vec(coverage(msa)))
```

It is possible to use this `DataFrame` and `Plots` to plot the sequence coverage of the MSA
and also an histogram of the number of sequences in each cluster:  

```@example msa_clusters
using StatsPlots # Plotting DataFrames
h = @df df histogram(:cluster, ylabel="nseq")
p = @df df plot(:cluster, :coverage, linetype=:scatter)
plot(p, h, nc=1, xlim=(0, nclusters(clusters)+1 ), legend=false)
png("msa_clusters_ii.png") # hide
nothing # hide
```  

![](msa_clusters_ii.png)  

We use the *Split-Apply-Combine* strategy, though the `by` function of the `DataFrames`
package, to select the sequence of highest coverage for each cluster.  

```@example msa_clusters
maxcoverage = by(df, :cluster, cl -> cl[ findmax(cl[:coverage])[2] ,
                 [:seqnum, :seqname, :coverage]])
```

```@example msa_clusters
p = @df maxcoverage plot(:cluster, :coverage, linetype=:scatter)
h = @df maxcoverage histogram(:cluster, ylabel="nseq")
plot(p, h, nc=1, xlim=(0, nclusters(clusters)+1 ), legend=false)
png("msa_clusters_iii.png") # hide
nothing # hide
```  

![](msa_clusters_iii.png)  

We can easily generate a mask using list comprehension, to select only the representative
sequences of the MSA (deleting the rest of the sequences with `filtersequences!`).  

```@example msa_clusters
cluster_references = Bool[ seqnum in maxcoverage[:seqnum] for seqnum in 1:nsequences(msa) ]
```

```@example msa_clusters
filtersequences!(msa, cluster_references)
```

```@example msa_clusters
plot(msa)
png("msa_clusters_iv.png") # hide
nothing # hide
```  

![](msa_clusters_iv.png)  
