```@meta
CurrentModule = MIToS.Information
```

```@setup log
@info "Information docs"
```

# [Information](@id Module-Information)

The `Information` module of MIToS defines types and functions useful to calculate
information measures (e.g. *Mutual Information* (MI) and *Entropy*) over a Multiple
Sequence Alignment (MSA). This module was designed to count `Residue`s
(defined in the `MSA` module) in special contingency tables (as fast as possible) and to
derive probabilities from these counts. Also, includes methods for applying corrections
to those tables, e.g. pseudocounts and pseudo frequencies. Finally, `Information` allows
to use these probabilities and counts to estimate information measures and other
frequency based values.  

```julia
using MIToS.Information # to load the Information module
```  

## Features

- Estimate multi dimensional frequencies and probability tables from sequences, MSAs, etc...
- Correction for small number of observations
- Correction for data redundancy on a MSA
- Estimate information measures
- Calculate corrected mutual information between residues  

## Contents

```@contents
Pages = ["Information.md"]
Depth = 4
```  

## Counting residues

MIToS Information module defines a multidimensional `ContingencyTable` type and two types
wrapping it, `Counts` and `Probabilities`, to store occurrences or probabilities.
The `ContingencyTable` type stores the contingency matrix, its marginal values and total.
These types are parametric, taking three ordered parameters:

- `T` : The type used for storing the counts or probabilities, e.g. `Float64`. It's
possible to use `BigFloat` if more precision it's needed.
- `N` : It's the dimension of the table and should be an `Int`.
- `A` : This should be a type, subtype of `ResidueAlphabet`, i.e.: `UngappedAlphabet`,
`GappedAlphabet` or `ReducedAlphabet`.

!!! note
    `ContingencyTable` can be used for storing probabilities or counts. The wrapper types
    `Probabilities` and `Counts` are mainly intended to dispatch in methods that need to
    know if the matrix has probabilities or counts, e.g. `entropy`. In general, the use of
    `ContingencyTable` is recommended over the use of `Probabilities` and `Counts`.

In this way, a matrix for storing pairwise probabilities of residues (without gaps) can be
initialized using:  

```@example inf_zeros
using MIToS.Information

Pij = ContingencyTable(Float64, Val{2}, UngappedAlphabet())
```  

**[High level interface]** It is possible to use the functions `count` and `probabilities`
to easily calculate the frequencies of sequences or columns of a MSA, where the number of
sequences/columns determine the dimension of the resulting table.  

```@example inf_count
using MIToS.Information
using MIToS.MSA # to use res"..." to create Vector{Residue}

column_i = res"AARANHDDRDC-"
column_j = res"-ARRNHADRAVY"
#   Nij[R,R] =   1     1   = 2

Nij = count(column_i, column_j)
```  

You can use `sum` to get the stored total:  

```@example inf_count
sum(Nij) # There are 12 Residues, but 2 are gaps
```  

Contingency tables can be indexed using `Int` or `Residue`s:  

```@example inf_count
Nij[2, 2] # Use Int to index the table
```  
```@example inf_count
Nij[Residue('R'), Residue('R')] # Use Residue to index the table
```  

!!! warning  
    The number makes reference to the specific index in the table e.g `[2,2]` references
    the second row and the second column. The use of the number used to encode the residue
    to index the table is dangerous. The equivalent index number of a residue depends on
    the used alphabet and `Int(Residue('X'))` will be always out of bounds.  

Indexing with `Residue`s works as expected. It uses the alphabet of the contingency table
to find the index of the `Residue`.

```@example inf_reduced
using MIToS.Information
using MIToS.MSA

alphabet = ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")

column_i = res"AARANHDDRDC-"
column_j = res"-ARRNHADRAVY"
#   Fij[R,R] =   1  1  1   = 3 # RHK

Fij = count(column_i, column_j, alphabet=alphabet)
```  
```@example inf_reduced
Fij[Residue('R'), Residue('R')] # Use Residue to index the table
```  

The function `getcontingencytable` allows to access the wrapped `ContingencyTable` in a
`Counts` object. You can use it, in combination with `normalize` to get a contingency table
of probabilities. The result can be wrapped inside a `Probabilities` object:  

```@example inf_reduced
Probabilities(normalize(getcontingencytable(Fij)))
```

#### Example: Plotting the probabilities of each residue in a sequence

Similar to the `count` function, the `probabilities` function can take at least one
sequence (vector of residues) and returns the probabilities of each residue. Optionally,
the keyword argument `alphabet` could be used to count some residues in the same cell
of the table.  

```@example inf_reduced
probabilities(res"AARANHDDRDC", alphabet=alphabet)
```

Here, we are going to use the `probabilities` function to get the residue probabilities of a
particular sequence from *UniProt*.

use the `getsequence` function, from the `MSA` module, to get the sequence from a `FASTA` downloaded from UniProt.  

```@repl
using MIToS.Information # to use the probabilities function
using MIToS.MSA # to use getsequence on the one sequence FASTA (canonical) from UniProt
seq = read("http://www.uniprot.org/uniprot/P29374.fasta", FASTA) # Small hack: read the single sequence as a MSA
probabilities(seq[1,:]) # Select the single sequence and calculate the probabilities
```  

!!! note
    In the previous example, using `getsequence(seq,1)` instead of `seq[1,:]` will return
    the sequence as a matrix with a single column to keep information for both dimensions.
    To use `probabilities` (or `count`) you can make use of the Julia's `vec` function to
    transform the matrix to a vector, e.g.: `probabilities(vec(getsequence(seq,1)))`.

```@setup inf_plotfreq
@info "Information: Plots"
using Plots
gr(size=(600,300))
using MIToS.Information # to use the probabilities function
using MIToS.MSA # to use getsequence on the one sequence FASTA (canonical) from UniProt
seq = read("http://www.uniprot.org/uniprot/P29374.fasta", FASTA) # Small hack: read the single sequence as a MSA
frequencies = probabilities(seq[1,:]) # Select the single sequence and calculate the probabilities
```  

```@example inf_plotfreq
using Plots # We choose Plots because it's intuitive, concise and backend independent
gr(size=(600,300))
```  

You can plot together with the probabilities of each residue in a given sequence, the
probabilities of each residue estimated with the BLOSUM62 substitution matrix. That matrix
is exported as a constant by the `Information` module as `BLOSUM62_Pi`.  

```@example inf_plotfreq
bar(
    1:20,
    [ frequencies  BLOSUM62_Pi ],
    lab = [ "Sequence"  "BLOSUM62"   ],
    alpha=0.5
    )
png("inf_plotfreq.png") # hide
nothing # hide
```  

![](inf_plotfreq.png)  

## Low count corrections

Low number of observations can lead to sparse contingency tables, that lead to wrong
probability estimations. It is shown in
[*Buslje et. al. 2009*![](./assets/external-link.png)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2672635/)
that low-count corrections, can lead to improvements in the contact prediction capabilities
of the Mutual Information. The Information module has available two low-count corrections:  

1. [Additive Smoothing![](./assets/external-link.png)](https://en.wikipedia.org/wiki/Additive_smoothing); the constant value pseudocount described in [*Buslje et. al. 2009*![](./assets/external-link.png)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2672635/).  
2. BLOSUM62 based pseudo frequencies of residues pairs, similar to [*Altschul et. al. 1997*![](./assets/external-link.png)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC146917/).  

```@example inf_msa
using MIToS.MSA

msa = read("http://pfam.xfam.org/family/PF09776/alignment/full", Stockholm)

filtercolumns!(msa, columngapfraction(msa) .< 0.5) # delete columns with 50% gaps or more

column_i = msa[:,1]
column_j = msa[:,2]
```

If you have a preallocated `ContingencyTable` you can use `count!` to fill it, this prevent
to create a new table as `count` do. However, you should note that `count!` **adds the new
counts to the pre existing values**, so in this case, we want to start with a table
initialized with zeros.  


```@example inf_msa
using MIToS.Information

const alphabet = ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")

Nij = ContingencyTable(Float64, Val{2}, alphabet)
```  

```@example inf_msa
#      table  weights         pseudocount      sequences...
count!(Nij,   NoClustering(), NoPseudocount(), column_i, column_j)
```  

!!! note
    You can use `NoClustering()` in places where clustering weights are required to not use
    weights. Also, `NoPseudocount()` in places where pseudocount values are required to not
    use pseudocounts.

In cases like the above, where there are few observations, it is possible to apply a
constant pseudocount to the counting table.  This module defines the type
`AdditiveSmoothing` and the correspond `fill!` and  `apply_pseudocount!` methods to
efficiently add or fill with a constant value each element of the table.

```@example inf_msa
apply_pseudocount!(Nij, AdditiveSmoothing(1.0))
```

**[High level interface.]** The `count` function has a `pseudocounts` keyword argument that
can take a `AdditiveSmoothing` value to easily calculate occurrences with pseudocounts. Also
the alphabet keyword argument can be used to chage the default alphabet (i.e. )


```@example inf_msa
count(column_i, column_j, pseudocounts=AdditiveSmoothing(1.0), alphabet=alphabet)
```  

To use the conditional probability matrix `BLOSUM62_Pij` in the calculation of pseudo
frequencies $G$ for the pair of residues $a$, $b$, it should be calculated first the real
frequencies/probabilities $p_{a,b}$. The observed probabilities are then used to estimate
the pseudo frequencies.  

$$G_{ab} = \sum_{cd}  p_{cd} \cdot BLOSUM62( a | c ) \cdot BLOSUM62( b | d )$$  

Finally, the probability $P$ of each pair of residues $a$, $b$ between the columns
$i$, $j$ is the weighted mean between the observed frequency $p$ and BLOSUM62-based
pseudo frequency $G$, where α is generally the number of clusters or the number of
sequences of the MSA and β is an empiric weight value. β was determined to be close
to `8.512`.  

$$P_{ab} = \frac{\alpha \cdot p_{ab} + \beta \cdot G_{ab} }{\alpha + \beta}$$

This could be easily achieved using the `pseudofrequencies` keyword argument of the
`probabilities` function. That argument can take a `BLOSUM_Pseudofrequencies` object that
is created with α and β as first and second argument, respectively.

```@example inf_msa
Pij = probabilities(column_i, column_j, pseudofrequencies=BLOSUM_Pseudofrequencies(nsequences(msa), 8.512))
```

You can also use `apply_pseudofrequencies!` in a previously filled probability contingency
table. i.e. `apply_pseudofrequencies!(Pij, BLOSUM_Pseudofrequencies(α, β))`

!!! warning
    `BLOSUM_Pseudofrequencies` can be only be applied in **normalized/probability** tables
    with `UngappedAlphabet`.  

## Correction for data redundancy in a MSA  

A simple way to reduce redundancy in a MSA without losing sequences, is clusterization and
sequence weighting. The weight of each sequence should be 1/N, where N is the number of
sequences in its cluster. The `Clusters` type of the `MSA` module stores the
weights. This vector of weights can be extracted (with the `getweight` function) and used
by the `count` and `probabilities` functions with the keyword argument `weights`. Also it's
possible to use the `Clusters` as second argument of the function `count!`.  


```@example inf_msa
clusters = hobohmI(msa, 62) # from MIToS.MSA
```

```@example inf_msa
count(msa[:,1], msa[:,2], weights=clusters)
```

## Estimating information measures on an MSA

The `Information` module has a number of functions defined to calculate information
measures from `Counts` and `Probabilities`:

- `entropy` : Shannon entropy (H)
- `marginal_entropy` : Shannon entropy (H) of the marginals
- `kullback_leibler` : Kullback-Leibler (KL) divergence
- `mutual_information` : Mutual Information (MI)
- `normalized_mutual_information` : Normalized Mutual Information (nMI) by Entropy
- `gap_intersection_percentage`
- `gap_union_percentage`

Information measure functions take optionally the base as the last positional argument
(default: `e`). You can use `2.0` to measure information in bits.

```@example inf_information
using MIToS.Information
using MIToS.MSA

Ni = count(res"PPCDPPPPPKDKKKKDDGPP") # Ni has the count table of residues in this low complexity sequence

H = entropy(Ni) # returns the Shannon entropy in nats (base e)
```

```@example inf_information
H = entropy(Ni, 2.0) # returns the Shannon entropy in bits (base 2)
```

Information module defines special iteration functions to easily and efficiently compute a
measure over a MSA. In particular, `mapcolfreq!` and `mapseqfreq!` map a function that takes
a table of `Counts` or `Probabilities`. The table is filled in place with the counts or
probabilities of each column or sequence of a MSA, respectively. `mapcolpairfreq!` and
`mapseqpairfreq!` are similar, but they fill the table using pairs of columns or sequences,
respectively.  

This functions take three positional arguments: the function `f` to be calculated, the
`msa` and `table` of `Counts` or `Probabilities`.  

After that, this function takes some keyword arguments:

- `weights` (default: `NoClustering()`) : Weights to be used for table counting.
- `pseudocounts` (default: `NoPseudocount()`) : `Pseudocount` object to be applied to table.
- `pseudofrequencies` (default: `NoPseudofrequencies()`) : `Pseudofrequencies` to be
applied to the normalized (probabilities) table.  

`mapcolpairfreq!` and `mapseqpairfreq!` also have a fourth positional argument `usediagonal`
that indicates if the function should be applied to identical element pairs
(default to `Val{true}`). This two functions also have an extra keyword argument
`diagonalvalue` (default to zero) to indicate the value used to fill the diagonal elements
if `usediagonal` is `Val{false}`.  

#### Example: Estimating *H(X)* and *H(X, Y)* over an MSA

In this example, we are going to use `mapcolfreq!` and `mapcolpairfreq!` to estimate
Shannon `entropy` of MSA columns *H(X)* and the joint entropy *H(X, Y)* of columns pairs,
respectively.  

```@setup inf_entropy
@info "Information: Entropy"
using Plots
gr()
```

```@example inf_entropy
using MIToS.MSA

msa = read("http://pfam.xfam.org/family/PF09776/alignment/full", Stockholm)
```

We are going to count residues to estimate the entropy. The `entropy` estimation is
performed over a rehused `Counts` object. The result will be a vector containing the
values estimated over each column without counting gaps (`UngappedAlphabet`).  

```@example inf_entropy
using MIToS.Information

Hx = mapcolfreq!(entropy, msa, Counts(ContingencyTable(Float64, Val{1}, UngappedAlphabet())))
```  

If we want the **joint entropy** between columns pairs, we need to use a bidimensional
table of `Counts` and `mapcolpairfreq!`.

```@example inf_entropy
Hxy = mapcolpairfreq!(entropy, msa, Counts(ContingencyTable(Float64, Val{2}, UngappedAlphabet())))
```  

In the above examples, we indicate the type of each occurrence in the counting and the probability table to use. Also, it's possible for some measures as **entropy** and **mutual information**, to estimate the values only with the count table (without calculate the probability table). Estimating measures only with a `ResidueCount` table, when this is possible, should be faster than using a probability table.  


```@example inf_entropy
Time_Pab = map(1:100) do x
    time = @elapsed mapcolpairfreq!(entropy, msa, Probabilities(ContingencyTable(Float64, Val{2}, UngappedAlphabet())))
end

Time_Nab = map(1:100) do x
    time = @elapsed mapcolpairfreq!(entropy, msa, Counts(ContingencyTable(Float64, Val{2}, UngappedAlphabet())))
end

using Plots
gr()

histogram( [Time_Pab Time_Nab],
    labels = ["Using ResidueProbability" "Using ResidueCount"],
    xlabel = "Execution time [seconds]" )

png("inf_entropy.png") # hide
nothing # hide
```  

![](inf_entropy.png)   

## Corrected Mutual Information  

MIToS ships with two methods to easily calculate corrected mutual information.  
The first is the algorithm described in [*Buslje et. al. 2009*![](./assets/external-link.png)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2672635/).
This algorithm can be accessed through the `buslje09` function and includes:  

1. Low count correction using `AdditiveSmoothing`
2. Sequence weighting after a `hobohmI` clustering
3. Average Product Correction (APC) proposed by
[Dunn et. al. 2008![](./assets/external-link.png)](http://bioinformatics.oxfordjournals.org/content/24/3/333),
through the `APC!` function that takes a MI matrix.
4. Z score correction using the functions `shuffle!` from the MSA module and `zscore`
from the `PairwiseListMatrices` package.  

```@docs
buslje09
```

The second, implemented in the `BLMI` function, has the same corrections that the above
algorithm, but use BLOSUM62 pseudo frequencies. This function is **slower** than
`buslje09` (at the same number of samples), but gives **better performance**
(for structural contact prediction) when the MSA has **less than 400 clusters** after a
Hobohm I at 62% identity.  

```@docs
BLMI
```

#### Example: Estimating corrected MI from an MSA

```@setup inf_buslje09
@info "Information: MI"
using Plots
gr()
```

```@example inf_buslje09
using MIToS.MSA
using MIToS.Information

msa = read("http://pfam.xfam.org/family/PF16078/alignment/full", Stockholm)
ZMIp, MIp  = buslje09(msa)
ZMIp
```

```@example inf_buslje09
ZBLMIp, BLMIp  = BLMI(msa)
ZBLMIp
```

## Visualize Mutual Information

You can use the function of the `Plots` package to visualize the Mutual Information (MI)
network between residues. As an example, we are going to visualize the MI between residues
of the Pfam domain *PF16078*. The `heatmap` is the simplest way to visualize the values of
the Mutual Information matrix.  

```@example inf_buslje09
using Plots
gr()

heatmap(ZMIp, yflip=true)
png("inf_heatmap.png") # hide
nothing # hide
```  

![](inf_heatmap.png)   

ZMIp is a Z score of the corrected MIp against its distribution on a random MSA
(shuffling the residues in each sequence), so pairs with highest values are more likely
to co-evolve. Here, we are going to use the top 1% pairs of MSA columns.  

```@example inf_buslje09
using PairwiseListMatrices # to use getlist
using Statistics # to use quantile

threshold = quantile(getlist(ZMIp), 0.99)
```

```@example inf_buslje09
ZMIp[ ZMIp .< threshold ] .= NaN
heatmap(ZMIp, yflip=true)
png("inf_heatmap_top.png") # hide
nothing # hide
```  

![](inf_heatmap_top.png)   

We are going to calculate the cMI (cumulative mutual information) value of each node.
Where cMI is a mutual information score per position that characterizes the extent of
mutual information "interactions" in its neighbourhood. This score is calculated as the
sum of MI values above a certain threshold for every amino acid pair where the particular
residue appears. This value defines to what degree a given amino acid takes part in a
mutual information network and we are going to indicate it using the node color.
To calculate cMI we are going to use the `cumulative` function:   

```@example inf_buslje09
cMI = cumulative(ZMIp, threshold)
```

```@setup comment_block
# # Setup block to hide this until PlotRecipes get fixed

# The nodes have an order, because they are columns in a MSA. So, the **arc diagram** it's
# useful to visualize long and short association between MSA positions. In general, long
# interactions has more interest.

# ` ` `@example inf_buslje09
# using PlotRecipes

# graphplot(ZMIp, size=(600,250), method=:arcdiagram) # , zcolor=cMI)
# png("inf_arcdiagram.png") # hide
# nothing # hide
# ` ` `  

# ![](inf_arcdiagram.png)   

# You can also use a **chord diagram** to see the same pattern.  

# ` ` `@example inf_buslje09
# graphplot(ZMIp, size=(600,600), method=:chorddiagram)
# png("inf_chorddiagram.png") # hide
# nothing # hide
# ` ` `  

# ![](inf_chorddiagram.png)   

```
