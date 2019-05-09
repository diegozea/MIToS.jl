## MIToS.jl Release Notes

### Changes from v2.3.0 to v2.4.0

MIToS v2.4 uses `Project.toml` and includes several bug fixes.

* The `SIFTS` module includes the `dbEnsembl` database and `warn`s again about unused databases.

### Changes from v2.2.0 to v2.3.0

MIToS v2.3 requires Julia v0.7 or v1.0. This release drops Julia 0.6 support.

* `Formatting.jl` is used in place of `Format.jl`.

* `SIFTS.get` returns the desired object or `missing` instead of `Nullable`s.

* `SIFTS` function doesn't `warn` about unused databases.

#### Julia 0.7/1.0 deprecations

* `bits` was deprecated to `bitstring`.

* `'` and `.'` are deprecated for alignments and sequences, use `transpose` or
`permutedims` instead. `ctranspose` is not longer available for matrices of `Residue`s.

### Changes from v2.1.2 to v2.2

* `PIR` `FileFormat` is included to read and write alignments in PIR/NBRF format.

* `Utils.Format` was renamed to `Utils.FileFormat`.

* `HTTP.jl` is used in place of `FTPClient.jl` and the deprecated `Requests.jl` in
`Utils.download_file` to download files.

* `Format.jl` is used in place of `Formatting.jl`.

* Solve bug in the printing of matrices of `Residue`s using `FileFormat`s.

### Changes from v2.1.1 to v2.1.2

* `FTPClient.jl` is used in `Utils.download_file` to download files from FTP.

* `CodecZlib.jl` is used in place of `GZip.jl` speeding up the parsing of compressed files.

* Improvements in MSA and PDB parsing speed.

* Improvement in `MSA.percentidentity` speed.

* `Information.gaussdca` now uses Julia's `serialize` and `deserialize` instead of `JLD`.

* `ROCAnalysis.jl` is not longer a dependency and it's now used with `@require` from
`Requires.jl`. To use the `AUC` function you need to do `using ROCAnalysis`.

### Changes from v2.1 to v2.1.1

* The script `Conservation.jl` was added to measure residue conservation of MSA columns.

* The script `SplitStockholm.jl` now has a progress bar thanks to Ellis Valentiner
@ellisvalentiner.

### Changes from v2.0 to v2.1

MIToS v2.1 requires Julia v0.6. This release drops Julia 0.5 support.

* `get_n_words(...` doesn't remove the last newline character, use `get_n_words(chomp(...`
to get the previous behaviour.

### Changes from v1.2.3 to v2.0

**MIToS 2.0** is the first MIToS version with **Julia 0.5** support
(It drops Julia 0.4 support). The last Julia version introduces new awesome features like
native multi-threading support, fast anonymous functions, generator expressions and more.
Also, the Julia package ecosystem has grown. So, MIToS was slightly redesigned to take
advantage of the new Julia capabilities. As a consequence, this version introduces several
breaking changes and new features.

##### Utils module

* `deleteitems!(vector::Vector, items)` is deprecated in favor of
`filter!(x -> x ∉ items, vector)`.

* `All` is used instead of MIToS 1.0 `"all"` or `"*"`, because it's possible to dispatch on it.

###### Vectorized queries are deprecated

Previous version of Utils included methods and types in order to overcome the performance
cost of functional programing in previous Julia versions. In particular, vectorized queries
were performed using subtypes of `AbstractTest`, in particular the `TestType`s `Is` and
`In` and the `TestOperation` `Not`. This types were used as argument to the query methods
`capture` and `isobject`. This operation were fused and vectorized with the methods:
`findobjects`, `collectobjects` and `collectcaptures`. All these functions and types are
deprecated in MIToS 2.0. Functional programming in Julia 0.5 is fast, so these methods
can be easily replace by Julia higher order functions like `find` and `filter` and lambda
expressions (anonymous functions).

##### MSA module

* `Residue` is now encoded as `Int` instead of being encoded as `UInt8`, allowing faster
indexation using `Int(res::Residue)`. More memory is used, since the residues are encoded
using 32 or 64 bits instead of 8 bits.

* `XAA` is now used to indicate unknown, ambiguous and non standard residues instead of `GAP`.

* Conversions to and from `UInt8` aren't supported now.

* More `Base` methods are extended to work with `Residue`: `bits`, `zero`, `one`
and `isvalid`.

* `empty(Annotations)` was deprecated, use `Annotations()` instead.

* `msa["seq_name",:]` now returns a `NamedArray{Residue,1}` instead of an aligned sequence,
use `getsequence(msa,"seqname")` to get an aligned sequence with annotations.

* The `names` function was replaced by the `sequencenames` function. A `columnnames`
function was also added.  

* Aligned sequences don't drop dimensions, so there are matrices instead of vectors. You can
use `vec(...)` or `squeeze(...,1)` to get a vector instead of the matrix.

* Indexing MSA objects with only one string is deprecated, use `msa["seqname",:]` instead
of `msa["seqname"]`.

* `empty!` doesn't take MSA objects anymore.

* `asciisequence` was replaced by `stringsequence`.

* `deletenotalphabetsequences` and the parse/read keyword argument `checkalphabet` are
deprecated since MIToS 2.0 uses Residue('X') to represent residues outside the alphabet. You
can use `filtersequences!(msa, vec(mapslices(seq -> !in(XAA, seq), msa, 2)))` to delete
sequences with unknown, ambiguous or non standard residues.

* `parse`/`read` and MSA file returns an `AnnotatedMultipleSequenceAlignment` by default.

* `shuffle_...columnwise!` and `shuffle_...sequencewise!` functions were deprecated in
favor of `shuffle!` and `shuffle` functions.

* `SequenceClusters` was renamed to `Clusters`.

* Residue alphabet types were added. All alphabet types are subtypes of `ResidueAlphabet`.
In particular, three types are exported: `GappedAlphabet`, `UngappedAlphabet` and
`ReducedAlphabet`. The last type allows the creation of custom reduced alphabets.

* In order to keep the sequence name, `AlignedSequence` and `AnnotatedAlignedSequence` are
now matrices instead of vectors.

##### PDB module

* The keyword argument `format` of `downloadpdb` should be a type (`PDBFile` or `PDBML`)
instead of a string (`pdb` or `xml`) as in MIToS 1.0.

* `read` and `parse` now has the `occupancyfilter` keyword argument.

* `read` and `parse` now has the `label` keyword argument for `PDBML` files.

* `residues`, `àtoms` and similiar functions don't take vectors or sets anymore. Use an
anonymous function instead, e.g.: `x -> x in set_of_residue_numbers`.

* The functions `isresidue`, `isatom` and `residuepairsmatrix` were added.

##### SIFTS module

* The `get` function has a more complex signature for `SIFTSResidue`s to make simpler
the access of data.

* `find`, `filter` and `filter` now takes a database type as a third parameter when a vector
of `SIFTSResidue`s is the second parameter. It allows to use a function that directly
operates over the database type if it's available.

* `SIFTSResidue`s now also store secondary structure data in the `sscode` and `ssname` fields.

##### Information module

* `ResidueProbability` and `ResidueCount` were deprecated in favor of `ContingencyTable`.
`Probabilities` and `Counts` were added as wrappers of `ContingencyTable` to allow dispach
in a some functions, e.g. `entropy`.

* The last parameter of contingency tables is now a subtype of `ResidueAlphabet` instead
of a `Bool`, i.e.: `UngappedAlphabet`, `GappedAlphabet` or `ReducedAlphabet`.

* Creation of empty contingecy tables chaged.
e.g. `zeros(ResidueProbability{Float64, 2, false})` changed to
`ContingencyTable(Float64, Val{2}, UngappedAlphabet())` and
`ResidueProbability{Float64, 2, false}()` changed to
`ContingencyTable{Float64, 2, UngappedAlphabet}(UngappedAlphabet())`.

* `count!` and `probabilities!` signatures changed. The first argument is alway a
`ContingencyTable`, the second positional argument a clustering weight object
(use `NoClustering()` to skip it), the third positional argument is a pseudocount object
(use `NoPseudocount()` to avoid the use of pseudocounts) and `probabilities!` takes also a
`Pseudofrequencies` object (use `NoPseudofrequencies()` to avoid pseudofrequencies). The
last positional arguments are the vector of residues used to fill the contingency table.

* `count` and `probabilities` now takes the sequences as only positional arguments. The
output is always a table of `Float64`. Both functions take the keyword arguments
`alphabet`, `weights` and `pseudocounts`. `probabilities` also has a `pseudofrequencies`
keyword argument.

* `apply_pseudofrequencies!` changed its signature. Now it takes a `ContingencyTable` and
a `Pseudofrequencies` object.

* The function `blosum_pseudofrequencies!` was deprecated in favor of introducing a
`BLOSUM_Pseudofrequencies` type as subtype of `Pseudofrequencies` to be used in
`probabilities`, `probabilities!` and `apply_pseudofrequencies!`.

* Because higher-order function are fast in Julia 0.5, measure types
(i.e. subtypes of `AbstractMeasure`) were deprecated in favor of functions. In particular,
`MutualInformation` was replaced with the `mutual_information` function,
`MutualInformationOverEntropy` was replaced with `normalized_mutual_information`,
`KullbackLeibler` was replaced with `kullback_leibler` and `Entropy` was replaced with
`entropy`.

* The functions `estimate`, `estimate_on_marginal` , `estimateincolumns` and
`estimateinsequences` were deprecated because measure types are not longer used.

* `estimate_on_marginal(Entropy...` was deprecated in favor of the `marginal_entropy`
function.

* `estimateincolumns` and `estimateinsequences` were deprecated in favor of `mapcolfreq!`,
`mapseqfreq!`, `mapcolpairfreq!` and `mapseqpairfreq`.

* Keyword argument `usegaps` is deprecated in `buslje09` and `BLMI` in favor of `alphabet`.

* `cumulative` function was added to calculate cumulative MI (cMI).

---

### Changes from v1.1 to v1.2.2

* `using Plots` to use `plot` with `AbstractVector{PDBResidue}` to visualize coordinates
of the C alpha of each residue.

* Re-exports `swap!` from **IndexedArrays.jl**.

* *[breaking change]* **Distances.jl** now uses `--inter` instead of `--intra`.

* *docs* and *cookbook* are now in [MIToSDocumentation](https://github.com/diegozea/MIToSDocumentation)

---

### Changes from v1.0 to v1.1

* **RecipesBase** is used to generate plot recipes for MIToS’ objects. MSA objects can be
visualized `using Plots` (thanks to Thomas Breloff @tbreloff ).

* Functions to perform structural superimposition were added to the `PDB` module
(thanks to Jorge Fernández de Cossío Díaz @cosio ) : `center!`, `kabsch`, `rmsd`.

* The `PDB` module adds the following functions to make easier structural comparison:
`getCA`, `CAmatrix`, `coordinatesmatrix`, `centeredcoordinates`, `centeredresidues`,
`change_coordinates`, `superimpose`, `mean_coordinates` and `rmsf`.

* When PDB or PDBML files are being parsed, It’s possible to indicate if only atoms with
the best occupancy should be loaded (`occupancyfilter=true`, `false` by default).

* When `PDBML` files are being parsed, is possible to used the new `label` keyword argument
to indicate if "auth" (`false`) or "label" (`true`) attributes should be used.

* `bestoccupancy!` was deprecated in favor of `bestoccupancy`.

* The `MSA` module export the function `percentsimilarity` to calculate the similarity
percent between aligned sequences.

* `msacolumn2pdbresidue` has two new keyword arguments, `strict` and `checkpdbname`, to
perform extra tests during the mapping between PDB and MSA residues.

* `msacolumn2pdbresidue` has a new `missings` keyword argument to indicate if missing
residues should be included in the mapping (default: `true`).

* The `MSA` now exports the `residue2three` and `three2residue` function to convert
`Residue`s to and from their three letter names.

* The `MSA` module now exports `sequencepairsmatrix`, `columnpairsmatrix`, `columnlabels`,
and `sequencelabels` to help in the construction of matrices for MSA sequences or columns
pairwise comparisons.

* The `Information` module, if `GaussDCA` is installed, allows to call its `gDCA` function
from MIToS through the `gaussdca` function.

* The `Information` module now exports the `KullbackLeibler` measure.

* Now is possible to `print` and `write` `PDBResidue`s as `PDBFile`s.

* The function `proximitymean` now has a keyword argument `include` to indicate if the
residue score should be included in the mean.

* The module `Scripts` inside the `Utils` module has a new function `readorparse` to help
parsing `STDIN` in MIToS’ scripts.

**MIToS v1.1** also includes several **bug fixes**, some **performance improvements** and a
more complete **documentation**.

---

### Changes from v0.1 to v1.0

* `Pfam` module for working with *Pfam* alignments and useful parameter optimization
functions (i.e. `AUC`).

* *[breaking change]* The `Clustering` module was deleted and its functions moved to the
`MSA` module.

* `MSA` uses `ClusteringResult` from the `Clustering.jl` package instead of `AbstractClusters`.

  * `Clusters` was renamed to `SequenceClusters`

  * `MSA` adds the `counts` and `assignments` functions from the `Clustering.jl` interface.

  * *[breaking change]* The `getnclusters` function is now `nclusters` in the `Clutering` module.

* *[breaking change]* All the MSA `...percentage` functions were renamed to `...fraction`
and `percent...` functions now return real percentages (not fractions) values.
Functions taking identity thresholds, now also take real percentages
(values between 0.0 and 100.0).

* *[breaking change]* Script command line arguments changed to: define the number of
workers, use STDIN and STDOUT (pipelines), get better output names, use real flag arguments.

* `InformationMeasure` renamed to `AbstractMeasure`.

* New functions added to `MSA` module.

  * `annotations`, `names`.

  * `meanpercentidentity` allows fast estimation of the mean percent identity between the sequences of a MSA.

* New function and type added to `Information` module.

  * `cumulative` to calculate cMI (cumulative mutual information) and similar cumulative scores.

  * `KullbackLeibler` to estimate conservation.

* `proximitymean` is defined in the `PDB` module to calculate pMI
(proximity mutual information) and other proximity scores.

* `contact` and `distance` have a vectorized form to create contact/distance maps.

* `NCol` file annotation with the number of columns in the original MSA.

* `BLMI` has `lambda` as a keyword argument for using additive smoothing.

* `BLMI` and `buslje09` accepts `samples=0` to avoid the Z score estimation.

* `read`/`parse` added the keyword argument `checkalphabet` for deleting sequences with non
standard amino acids.

* `read`/`parse` added the keyword argument `keepinserts` for keep insert columns
(It creates an `Aligned` column annotation).

**MIToS v1.0** also includes several **bug fixes** and a more complete **documentation**.
