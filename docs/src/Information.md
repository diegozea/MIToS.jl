```@meta
CurrentModule = MIToS.Information
```

```@setup log
@info "Information docs"
```

# [Information](@id Module-Information)

Extracting evolutionary signals, such as conservation and coevolution, from
Multiple Sequence Alignments (MSAs) is a common task in bioinformatics. There are
several methods to estimate these signals, including information measures like
*Shannon Entropy*—to assess the conservation of a position—and *Mutual Information*—to
assess the coevolution between two positions. The `Information` module of MIToS defines
types and functions useful for calculating those information measures over an MSA.
This module was designed to count `Residue`s (defined in the `MSA` module) in special
contingency tables (as fast as possible) and to derive probabilities from these counts.
It also includes methods for applying corrections to those tables, e.g., pseudo counts and
pseudo frequencies. Finally, `Information` allows using probabilities and counts
to estimate information measures and other frequency-based values.

```julia
using MIToS.Information # to load the Information module
```

## Features

  - Estimate multi-dimensional frequencies (counts) and probability tables from sequences,
    MSA columns, etc...
  - Corrections for a small number of observations
  - Corrections for data redundancy on an MSA
  - Estimate information measures such as Shannon entropy, mutual information, etc...
  - Calculate corrected mutual information between residues

## Contents

```@contents
Pages = ["Information.md"]
Depth = 4
```

## Counting residues

You can use the `frequencies` and `probabilities` functions to easily calculate the
amino acid frequencies or probabilities of a sequence or a column of an MSA. The number of
sequences/columns determines the dimension of the resulting table. Let's see an example
using the `frequencies` function to calculate the frequencies of each pair of residues in
two columns of an MSA. The resulting `Frequencies` object contains the bidimensional 
contingency table, the marginal values, and the total.

```@example inf_count
using MIToS.Information
using MIToS.MSA # to use res"..." to create Vector{Residue}

column_i = res"AARANHDDRDC-"
column_j = res"-ARRNHADRAVY"
#   Nij[R,R] =   1     1   = 2

Nij = frequencies(column_i, column_j)
```

You can use `sum` to get the stored total:

```@example inf_count
sum(Nij)
```

Here, the total is `10.0`, because there are `12` residues on each column, but `2` are gaps. Since the default alphabet used by `frequency` is `UngappedAlphabet()`, the gaps are 
not counted.

Contingency tables can be indexed using `Int` (the coordinates in the table)
or `Residue`s. For example, to get the number of times an arginine (*R*) is
found in the same sequence in both columns, you can use:

```@example inf_count
Nij[Residue('R'), Residue('R')]
```

Or, since the arginine is in the second row and the second column of the table—this is 
because the arginine is the second residue in the `UngappedAlphabet()`—you can use:

```@example inf_count
Nij[2, 2]
```

!!! warning "Indexing with `Int`"
    
    The index refers to the specific position in the table, e.g., `[2,2]` references the 
    second row and the second column. For `GappedAlphabet()` and `UngappedAlphabet()`, 
    the index is the same as the position of the residue in the alphabet. 
    However, for `ReducedAlphabet()`, the index will be the number of the group to which 
    the residue belongs. For example, if the alphabet is 
    `ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")`, the index of the arginine (*R*) 
    will be `3` rather than `2`. Therefore, **it is recommended** that you use a `Residue` 
    object to index the table to avoid subtle bugs if the alphabet changes.

You can change the alphabet used to count residues by setting the `alphabet` keyword
of the `frequencies` function. Let's see an example with an reduced alphabet:

```@example inf_reduced
using MIToS.Information
using MIToS.MSA

alphabet = ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")

column_i = res"AARANHDDRDC-"
column_j = res"-ARRNHADRAVY"
#   Fij[R,R] =   1  1  1   = 3 # RHK

Fij = frequencies(column_i, column_j, alphabet = alphabet)
```

```@example inf_reduced
Fij[Residue('R'), Residue('R')]
```

### Example: Plotting the probabilities of each residue in a sequence

Like the `frequencies` function, the `probabilities` function can take at least one
sequence or column (a vector of `Residue` objects) and return the probabilities of each 
residue. Optionally, as before, the keyword argument `alphabet` could be used to count 
some residues in the same table cell.

```@example inf_reduced
probabilities(res"AARANHDDRDC", alphabet = alphabet)
```

Here, we are going to use the `probabilities` function to get the residue probabilities of a
particular sequence from *UniProt*.

```@setup inf_plotfreq
@info "Information: Plots"
using Plots
gr(size=(600,300))
```

```@repl inf_plotfreq
using MIToS.Information
using MIToS.MSA
file_name = "http://www.uniprot.org/uniprot/P29374.fasta"
sequences = read_file(file_name, FASTASequences)
Pa = probabilities(sequences[1])
```

```@example inf_plotfreq
using Plots # We choose Plots because it's intuitive and concise
gr(size = (600, 300))
```

You can plot together with the probabilities of each residue in a given sequence, the
probabilities of each residue estimated with the BLOSUM62 substitution matrix. That matrix
is exported as a constant by the `Information` module as `BLOSUM62_Pi`.

```@example inf_plotfreq
bar(1:20, [Pa BLOSUM62_Pi], lab = ["Sequence" "BLOSUM62"], alpha = 0.5)
png("inf_plotfreq.png") # hide
nothing # hide
```

![](inf_plotfreq.png)

### Low-level interface

You can work directly with the `ContingencyTable` type if you want performance.
The MIToS Information module also defines two types wrapping this multi-dimensional array;
`Frequencies` and `Probabilities`, to store occurrences or probabilities, respectively.
The `ContingencyTable` type stores the contingency matrix, its marginal values, and its total.
These three types are parametric, taking three ordered parameters:

  - `T`: The type used for storing the counts or probabilities, e.g., `Float64`. If more  precision is needed, `BigFloat` is possible.
  - `N`: It's the dimension of the table and should be an `Int`.
  - `A`: This should be a type, subtype of `ResidueAlphabet`, i.e.: `UngappedAlphabet`, `GappedAlphabet`, or `ReducedAlphabet`.

!!! note

    `ContingencyTable` can be used to store probabilities or counts. The wrapper types
    `Probabilities` and `Frequencies` are mainly intended to dispatch in methods that need
    to know if the matrix has probabilities or counts. For example, the implementation of `shannon_entropy` is faster when the table is a `Frequencies` object.

In this way, a matrix for storing pairwise probabilities of residues (without gaps) can be
initialized using:

```@example inf_zeros
using MIToS.Information

Pij = ContingencyTable(Float64, Val{2}, UngappedAlphabet())
```

Note that the dimension of the table is indicated using `Val{2}` rather than `2`, so that Julia knows that the table is two-dimensional at compile time.

You can also obtain the `ContingencyTable` wrapped in a `Frequencies` object using 
the `getcontingencytable` function. For example, you can take a table of counts and 
normalize it using the `normalize` function to get a table of probabilities that you can 
later wrap in a `Probabilities` object:

```@example inf_convert
using MIToS.MSA
using MIToS.Information
column_i = res"AARANHDDRDC-"
column_j = res"-ARRNHADRAVY"
Fij = frequencies(column_i, column_j)
Pij = Probabilities(normalize(getcontingencytable(Fij)))
```

## Low count corrections

A low number of observations can lead to sparse contingency tables that lead to wrong
probability estimations. It is shown in [10.1093/bioinformatics/btp135](@citet)
that low-count corrections, can lead to improvements in the contact prediction capabilities
of the mutual information. The Information module has available two low-count corrections:

  1. [`AdditiveSmoothing`](@ref): [Additive smoothing![]
     (./assets/external-link.png)](https://en.wikipedia.org/wiki/Additive_smoothing) 
     corrects the frequencies by adding a constant value to each cell of the table.
     It can be used to represent the constant value pseudocount described in [10.1093/bioinformatics/btp135](@citet). For example: `AdditiveSmoothing(1.0)` can be used to 
     add `1.0` to each cell of the table.
  2. [`BLOSUM_Pseudofrequencies`](@ref): BLOSUM62-based pseudo frequencies for residues pairs, 
     similar to [10.1093/nar/25.17.3389](@citet) and described in 
     [10.1093/bioinformatics/btw646](@cite)

The `frequencies` and `probabilities` functions have a `pseudocounts` keyword argument that
can take an `AdditiveSmoothing` value to calculate occurrences and probabilities with 
pseudo counts.

```@example inf_pseudocounts
using MIToS.MSA
using MIToS.Information
seq = res"AAAVRNHDDRDCPPPGGPPG"
frequencies(seq, pseudocounts = AdditiveSmoothing(1.0))
```

The `probabilities` function can also take a `pseudofrequencies` keyword argument that can
take a `BLOSUM_Pseudofrequencies` object to calculate probabilities with pseudo 
frequencies. Note that these pseudofrequencies can only be used on bidimensional tables, 
i.e. passing only two sequences or columns to the `probabilities` function, and using `UngappedAlphabet()` (the default alphabet).

```@example inf_pseudofrequencies
using MIToS.MSA
using MIToS.Information
column_i = res"ARANHDDRDC"
column_j = res"-RRNHADRAV"
probabilities(column_i, column_j, pseudofrequencies = BLOSUM_Pseudofrequencies(10, 8.512))
```

### Low-level interface

If you have a preallocated `ContingencyTable`, for performance
reasons, you can use `frequencies!` or the `probabilities!` functions to fill it. That 
prevents the creation of a new table as `frequencies` and `probabilities` do. However, you
should note that `frequencies!` and `probabilities!` **add the new counts to the pre-existing values**. Therefore, you can use the `cleanup!` function to set the table to zero 
before filling it. 

Let's see an example, in this case we start with a table of zeros, `Nij`, and fill it 
with the counts of the residues in the columns `column_i` and `column_j` using the
`frequencies!` function. Since we start with a new table, we don't need to use `cleanup!`.

```@example inf_msa
using MIToS.MSA
using MIToS.Information

file_name = "https://raw.githubusercontent.com/diegozea/MIToS.jl/master/docs/data/PF18883.stockholm.gz"

msa = read_file(file_name, Stockholm)

column_i = msa[:, 1]
column_j = msa[:, 2]

const ALPHABET = ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")

Nij = ContingencyTable(Float64, Val{2}, ALPHABET)

frequencies!(Nij, column_i, column_j)
```

In cases like the above, where there are few observations, applying a constant pseudocount 
to the contingency table could be beneficial. This module defines the type
`AdditiveSmoothing` and the corresponding `fill!` and  [`apply_pseudocount!`](@ref) methods 
to efficiently fill or add a constant value to each element of the table.

```@example inf_msa
apply_pseudocount!(Nij, AdditiveSmoothing(1.0))
```


If we work with a normalized contingency table or a `Probabilities` object and use the 
`UnappedAlphabet()`, we can apply the BLOSUM62-based pseudo frequencies to the table. 
For that, we can use the `BLOSUM_Pseudofrequencies` type and the `apply_pseudofrequencies!` 
function. That function first needs to calculate the actual probabilities $p_{a,b}$ for 
each pair of residues $a$, $b$. Then, it uses the conditional probability matrix 
`BLOSUM62_Pij` and the observed probabilities to calculate pseudo frequencies $G$.

$$G_{ab} = \sum_{cd}  p_{cd} \cdot BLOSUM62( a | c ) \cdot BLOSUM62( b | d )$$

The `apply_pseudofrequencies!` function calculates the probability $P$ of each pair of 
residues $a$, $b$ between the columns $i$, $j$ as the weighted mean between the observed 
frequency $p$ and BLOSUM62-based pseudo frequency $G$. The former has α as weight, and the 
latter has the parameter β. We use the number of sequences or sequence clusters as α, and 
β is an empiric weight value that we determined to be close to `8.512`.

$$P_{ab} = \frac{\alpha \cdot p_{ab} + \beta \cdot G_{ab} }{\alpha + \beta}$$

The `BLOSUM_Pseudofrequencies` type is defined with two parameters, α and β, that are 
set using the positional arguments of the constructor.

```@example inf_pseudofrequencies
using MIToS.MSA
using MIToS.Information
column_i = res"ARANHDDRDC"
column_j = res"-RRNHADRAV"
Pij = ContingencyTable(Float64, Val{2}, UngappedAlphabet())
probabilities!(Pij, column_i, column_j)
α = 10
β = 8.512
apply_pseudofrequencies!(Pij, BLOSUM_Pseudofrequencies(α, β))
```

## Correction for data redundancy in a MSA

A simple way to reduce redundancy in an MSA without losing sequences is through 
clusterization and sequence weighting. The weight of each sequence should be $1/N$, 
where $N$ is the number of sequences in its cluster. The `Clusters` type of the `MSA` 
module stores the weights. This vector of weights can be extracted 
(with the `getweight` function) and used by the `frequencies` and `probabilities` functions 
with the keyword argument `weights`. Also, using the `Clusters` as the second argument of 
the function `frequencies!` is possible.

```@example inf_msa
clusters = hobohmI(msa, 62) # from MIToS.MSA
```

```@example inf_msa
frequencies(msa[:, 1], msa[:, 2], weights = clusters)
```

## Estimating information measures on an MSA

The `Information` module has a number of functions defined to calculate information
measures from `Frequencies` and `Probabilities`:

  - `shannon_entropy` : Shannon entropy (H)
  - `marginal_entropy` : Shannon entropy (H) of the marginals
  - `kullback_leibler` : Kullback-Leibler (KL) divergence
  - `mutual_information` : Mutual Information (MI)
  - `normalized_mutual_information` : Normalized Mutual Information (nMI) by Entropy
  - `gap_intersection_percentage`
  - `gap_union_percentage`

Information measure functions take optionally the base as a keyword argument (default: `ℯ`).
You can set `base=2` to measure information in bits.

```@example inf_information
using MIToS.Information
using MIToS.MSA

Ni = frequencies(res"PPCDPPPPPKDKKKKDDGPP") # Ni has the count table of residues in this low complexity sequence

H = shannon_entropy(Ni) # returns the Shannon entropy in nats (base e)
```

```@example inf_information
H = shannon_entropy(Ni, base = 2) # returns the Shannon entropy in bits (base 2)
```

Information module defines special iteration functions to easily and efficiently compute a
measure over a MSA. In particular, `mapcolfreq!` and `mapseqfreq!` map a function that takes
a table of `Frequencies` or `Probabilities`. The table is filled in place with the counts or
probabilities of each column or sequence of a MSA, respectively. `mapcolpairfreq!` and
`mapseqpairfreq!` are similar, but they fill the table using pairs of columns or sequences,
respectively.

This functions take three positional arguments: the function `f` to be calculated, the
`msa` and `table` of `Frequencies` or `Probabilities`.

After that, this function takes some keyword arguments:

  - `weights` (default: `NoClustering()`) : Weights to be used for table counting.
  - `pseudocounts` (default: `NoPseudocount()`) : `Pseudocount` object to be applied to table.
  - `pseudofrequencies` (default: `NoPseudofrequencies()`) : `Pseudofrequencies` to be applied to the normalized (probabilities) table.
  - `usediagonal` (default: `true`) : Indicates if the function should be applied to pairs containing the same sequence or column.
  - `diagonalvalue` (default to zero) : The value that fills the diagonal elements of the table if `usediagonal` is `false`.

#### Example: Estimating *H(X)* and *H(X, Y)* over an MSA

In this example, we are going to use `mapcolfreq!` and `mapcolpairfreq!` to estimate
Shannon `shannon_entropy` of MSA columns *H(X)* and the joint entropy *H(X, Y)* of columns pairs,
respectively.

```@setup inf_entropy
@info "Information: Entropy"
using Plots
gr()
```

```@example inf_entropy
using MIToS.MSA

msa = read_file(
    "https://raw.githubusercontent.com/diegozea/MIToS.jl/master/docs/data/PF18883.stockholm.gz",
    Stockholm,
)
```

We are going to count residues to estimate the Shannon entropy. The `shannon_entropy`
estimation is performed over a rehused `Frequencies` object. The result will be a vector
containing the values estimated over each column without counting gaps (`UngappedAlphabet`).

```@example inf_entropy
using MIToS.Information

Hx = mapcolfreq!(
    shannon_entropy,
    msa,
    Frequencies(ContingencyTable(Float64, Val{1}, UngappedAlphabet())),
)
```

If we want the **joint entropy** between columns pairs, we need to use a bidimensional
table of `Frequencies` and `mapcolpairfreq!`.

```@example inf_entropy
Hxy = mapcolpairfreq!(
    shannon_entropy,
    msa,
    Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet())),
)
```

In the above examples, we indicate the type of each occurrence in the counting and the probability table to use. Also, it's possible for some measures as **entropy** and **mutual information**, to estimate the values only with the count table (without calculate the probability table). Estimating measures only with a `ResidueCount` table, when this is possible, should be faster than using a probability table.

```@example inf_entropy
Time_Pab = map(1:100) do x
    time = @elapsed mapcolpairfreq!(
        shannon_entropy,
        msa,
        Probabilities(ContingencyTable(Float64, Val{2}, UngappedAlphabet())),
    )
end

Time_Nab = map(1:100) do x
    time = @elapsed mapcolpairfreq!(
        shannon_entropy,
        msa,
        Frequencies(ContingencyTable(Float64, Val{2}, UngappedAlphabet())),
    )
end

using Plots
gr()

histogram(
    [Time_Pab Time_Nab],
    labels = ["Using ResidueProbability" "Using ResidueCount"],
    xlabel = "Execution time [seconds]",
)

png("inf_entropy.png") # hide
nothing # hide
```

![](inf_entropy.png)

## Corrected Mutual Information

MIToS ships with two methods to easily calculate corrected mutual information.
The first is the algorithm described in [10.1093/bioinformatics/btp135](@citet).
This algorithm can be accessed through the `buslje09` function and includes:

 1. Low count correction using `AdditiveSmoothing`
 2. Sequence weighting after a `hobohmI` clustering [10.1002/pro.5560010313](@cite)
 3. Average Product Correction (APC) proposed by [10.1093/bioinformatics/btm604](@citet),
    through the `APC!` function that takes a MI matrix.
 4. Z score correction using the functions `shuffle_msa!` from the MSA module and `zscore`
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

msa = read_file(
    "https://raw.githubusercontent.com/diegozea/MIToS.jl/master/docs/data/PF18883.stockholm.gz",
    Stockholm,
)
ZMIp, MIp = buslje09(msa)
ZMIp
```

```@example inf_buslje09
ZBLMIp, BLMIp = BLMI(msa)
ZBLMIp
```

## Visualize Mutual Information

You can use the function of the `Plots` package to visualize the Mutual Information (MI)
network between residues. As an example, we are going to visualize the MI between residues
of the Pfam domain *PF18883*. The `heatmap` is the simplest way to visualize the values of
the Mutual Information matrix.

```@example inf_buslje09
using Plots
gr()

heatmap(ZMIp, yflip = true)
png("inf_heatmap.png") # hide
nothing # hide
```

![](inf_heatmap.png)

ZMIp is a Z score of the corrected MIp against its distribution on a random MSA
(shuffling the residues in each sequence), so pairs with highest values are more likely
to coevolve. Here, we are going to use the top 1% pairs of MSA columns.

```@example inf_buslje09
using PairwiseListMatrices # to use getlist
using Statistics # to use quantile

threshold = quantile(getlist(ZMIp), 0.99)
```

```@example inf_buslje09
ZMIp[ZMIp.<threshold] .= NaN
heatmap(ZMIp, yflip = true)
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
