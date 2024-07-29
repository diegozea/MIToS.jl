## MIToS.jl Release Notes

### Changes from v2.22.0 to v3.0.0

**MIToS v3.0.0** requires Julia v1.9 or higher, dropping support for older versions. This
release introduces several breaking changes to improve the usability of the package.
When possible, deprecation warnings are used to inform you of the changes.

#### MIToS.MSA

The MSA module now includes ways to read, write, and work with unaligned protein sequences:

  - The `MSA` module now exports the `AnnotatedSequence` type to represent a single protein
    sequence with annotations. This type is a subtype of the new `AbstractSequence` type,
    a subtype of the new `AbstractResidueMatrix` type.

  - The `MSA` module now exports the `sequence_id` function to get the identifier of a
    sequence object.
  - The `MSA` module now defines the `FASTASequences`, `PIRSequences`, and `RawSequences`
    file formats to read and write (unaligned) protein sequences in FASTA, PIR, and raw
    formats, respectively.
  - *[Breaking change]* The behavior of the `getannotresidue`, `getannotsequence`,
      `setannotresidue!`, and `setannotsequence!` functions have changed for sequences objects,
    such as `AnnotatedSequence`, `AnnotatedAlignedSequence`, and `AlignedSequence`. Now, these
    functions take the feature name, rather than the sequence name, as the second
    positional argument. As an example of migration,
      `getannotsequence(sequence, "sequence_name", "feature_name")` should be replaced by
      `getannotsequence(sequence, "feature_name")`. You still need to specify the sequence name
    when working with MSA objects.

Other changes in the MSA module are:

  - *[Breaking change]* The `join` function for `AnnotatedMultipleSequenceAlignment` objects
    is deprecated in favor of the `join_msas` function.

  - *[Breaking change]* The `Clusters` type is no longer a subtype of `ClusteringResult` from
    the `Clustering.jl` package. Instead, the `Clusters` type is now a subtype of the new
    `AbstractCluster` type. Support for the `Clustering.jl` interface is still available
    through package extensions. You now need to load the `Clustering.jl` package to use the
    `assignments`, `nclusters`, and `counts` functions.

#### MIToS.PDB

The PDB module now depends on the `BioStructures` package. The main changes in the PDB
module are:

  - The `PDB` module now exports the `MMCIFFile` file format to read and write PDB files in
    the mmCIF format (using `BioStructures` under the hood).

  - *[Breaking change]* The `download_alphafold_structure` function can now download the
    predicted structures from the *AlphaFold Protein Structure Database* using the mmCIF
    format (`format=MMCIFFile`). This is the new default format. Therefore, you should use
    `format=PDBFile` to get a PDB file like before. For example,
      `download_alphafold_structure("P00520")` in previous versions is the same as
      `download_alphafold_structure("P00520", format=PDBFile)` in this version.
  - *[Breaking change]* The `downloadpdb` function now returns a mmCIF file by default.
    Therefore, you should use `format=PDBML` to get a PDBML file. As an example of migration,
    `downloadpdb("1IVO")` should be replaced by `downloadpdb("1IVO", format=PDBML)`, unless
    you want to get a mmCIF file.
  - *[Breaking change]* The `PDBAtom` type now adds two extra fields: `alt_id` and `charge`
    to represent the alternative location indicator and the atom's charge, respectively.
    This improves the compatibility with the mmCIF format and the `BioStructures` package.
  - *[Breaking change]* The `query_alphafolddb` function now returns the EntrySummary object
    of the returned JSON response instead of the Root list. Therefore, there is no need to
    take the first element of the list to get the required information. For example,
      `query_alphafolddb("P00520")[1]["uniprotId"]` would be replaced by
      `query_alphafolddb("P00520")["uniprotId"]`.

#### MIToS.Utils.Scripts

  - *[Breaking change]* The `MIToS.Utils.Scripts` module and the MIToS scripts have been
    moved to their package at [MIToS_Scripts.jl](https://github.com/MIToSOrg/MIToS_Scripts.jl).
    Therefore, the `MIToS.Utils.Scripts` module is no longer exported. This allows for a
    reduction in the number of MIToS dependencies and improved load time.

### Changes from v2.21.0 to v2.22.0

This versions introduces several breaking changes to improve the usability of the
`Information` module. The main changes are:

  - *[Breaking change]* The `Information` module deprecates the `Counts` type in favor of
    the new `Frequencies` type. The new type as the same signature and behavior as the old one.

  - *[Breaking change]* The `count` function on sequences has been deprecated in favor of the
    `frequencies` function, which has the same signature and behavior as the old one.
  - *[Breaking change]* The `count!` function is deprecated in favor of `frequencies!`.
    The new function use keyword arguments to define the weights and pseudocounts. As an
    example of migration, `count!(table, weights, pseudocounts, seqs...)` should be replaced
    by `frequencies!(table, seqs..., weights=weights, pseudocounts=pseudocounts)`.
  - *[Breaking change]* The `probabilities!` method using positional arguments for the
    weights, pseudocounts and pseudofrequencies is deprecated in favor the one that uses
    keyword arguments. As an example of migration,
    `probabilities!(table, weights, pseudocounts, pseudofrequencies, seqs...)`
    should be replaced by
    `probabilities!(table, seqs..., weights=weights, pseudocounts=pseudocounts, pseudofrequencies=pseudofrequencies)`.
  - *[Breaking change]* The `Information` has deprecated the `entropy` method on
    `Frequencies` and `Probabilities` in favor of the `shannon_entropy` function. The
    definition of the base is now done using the `base` keyword argument. As an example of
    migration, `entropy(p, 2)` should be replaced by `shannon_entropy(p, base=2)`.
  - *[Breaking change]* The `marginal_entropy` methods based on positional arguments are
    deprecated in favor of a method relying on the `margin` and `base` keyword arguments.
    As an example of migration, `marginal_entropy(p, 2, 2.0)` should be replaced by
    `marginal_entropy(p, margin=2, base=2.0)`.
  - *[Breaking change]* The `mutual_information` method based on positional arguments is
    deprecated in favor of a method relying on the `base` keyword argument. As an example of
    migration, `mutual_information(p, 2)` should be replaced by `mutual_information(p, base=2)`.
  - *[Breaking change]* The `mapcolpairfreq!` and `mapseqpairfreq!` functions now uses the
    boolean `usediagonal` keyword argument to indicate if the function should be applied to
    the diagonal elements of the matrix (the default is `true`). Before, this was done passing
    `Val{true}` or `Val{false}` as the last positional argument.
  - The `mapcolfreq!`, `mapseqfreq!`, `mapcolpairfreq!`, and `mapseqpairfreq!` methods using
    keyword arguments, now pass the extra keyword arguments to the mapped function.
  - The `Information` module now exports the `mapfreq` function that offers a more high-level
    interface to the `mapcolfreq!`, `mapseqfreq!`, `mapcolpairfreq!`, and `mapseqpairfreq!`
    functions. This function allows the user to map a function to the residue frequencies or
    probabilities of the columns or sequences of an MSA. When `rank = 2`, the function is
    applied to pairs of sequences or columns.
  - The `Information` module now exports methods of the `shannon_entropy`, `kullback_leibler`,
    `mutual_information`, and `normalized_mutual_information` functions that take an
    `AbstractArray{Residue}` as input, e.g. an MSA. Those methods use the `mapfreq` function
    under the hood to ease the calculation of the information measures on MSAs.
  - The `frequencies!`, `frequencies`, `probabilities!`, and `probabilities` functions now
    accept arrays of `Residue`s of any dimension. Therefore, there is no need to use the
    `vec` function to convert the arrays to vectors.
  - The `MSA` module now exports the `WeightType` union type to represent `weights`.

### Changes from v2.20.0 to v2.21.0

  - *[Breaking change]* The `buslje09` and `BLMI` functions from the `Information` module does
    not longer accept a filename and a file format as arguments. You should explicitly read
    the MSA using the `read_file` function and then run the `buslje09` or `BLMI` functions
    on the returned MSA object. As an example of migration, `buslje09("msa.sto", "Stockholm")`
    should be replaced by `buslje09(read_file("msa.sto", Stockholm))`.

### Changes from v2.19.0 to v2.20.0

  - *[Breaking change]* The PDB module has deprecated `residues` and `@residues` in favor of
    the `select_residues` function that uses keyword arguments.
    So, `residues(pdb, "1", "A", "ATOM", All)` or `@residues pdb "1" "A" "ATOM" All` should be
    replaced by `select_residues(pdb, model="1", chain="A", group="ATOM")`.

  - *[Breaking change]* The PDB module has deprecated `atoms` and `@atoms` in favor of
    the `select_atoms` function that uses keyword arguments.
    So, `atoms(pdb, "1", "A", "ATOM", All, "CA")` or `@atoms pdb "1" "A" "ATOM" All "CA"` should be
    replaced by `select_atoms(pdb, model="1", chain="A", group="ATOM", atom="CA")`.
  - *[Breaking change]* The PDB module has deprecated the methods of the `isresidue` and
    `residuesdict` functions that rely on positional arguments in favor of the keyword arguments.
    So, `isresidue(pdb, "1", "A", "ATOM", "10")` should be replaced by
    `isresidue(pdb, model="1", chain="A", group="ATOM", residue="10")`. Similarly,
    `residuesdict(pdb, "1", "A", "ATOM", All)` should be replaced by
    `residuesdict(pdb, model="1", chain="A", group="ATOM")`.

### Changes from v2.18.0 to v2.19.0

  - *[Breaking change]* The `shuffle` and `shuffle!` functions are deprecated in favor of the
    `shuffle_msa` and `shuffle_msa!` functions. The new functions take `dims` and
    `fixedgaps` as keyword arguments instead of taking them as positional ones. The new
    functions add a last positional argument to allow the selection of specific sequences
    or columns to shuffle. Also, it adds the `fixed_reference` keyword argument to keep the
    residues in the reference sequence fixed during the shuffling. As an example of migration,
    `shuffle!(msa, 1, false)` should be replaced by `shuffle_msa!(msa, dims=1, fixedgaps=false)`.

### Changes from v2.17.0 to v2.18.0

  - *[Breaking change]* The `read`, `parse`, `write`, and `print` functions for different
    `FileFormat`s have been deprecated in favor of the `read_file`, `parse_file`,
    `write_file`, and `print_file` functions. The new functions keep the same signature and
    behavior as the old ones.

### Changes from v2.16.0 to v2.17.0

  - *[Breaking change]* The `download_file` now uses the `Downloads.jl` module instead of
    `HTTP.jl`. Therefore, the `download_file` function now accepts the `Downloads.download`
    keyword arguments. In particular, the `redirect` and `proxy` keyword arguments are no
    longer needed.

  - The `MSA` module now exports the `A2M` and `A3M` file formats, to allow reading and
    writing MSA files in these formats.

### Changes from v2.15.0 to v2.16.0

MIToS v2.16.0 drops support for *Julia 1.0*. This release requires *Julia 1.6* or higher.

  - *[Breaking change]* The `transpose` function is now deprecated for MSA and sequences
    (`AbstractAlignedObject`s). Use `permutedims` instead.

  - *[Breaking change]* MIToS is now using `JSON3.jl` instead of `JSON.jl`. That change the
    returned type of `getpdbdescription` from `Dict{String, Any}` to `JSON3.Object`.
    Since the `JSON3.Object` supports the `Dict` interface, the change should not cause any
    issues. If you want to convert the returned `JSON3.Object` to a `Dict{String, Any}`
    you can use the `MIToS.PDB.JSON3.copy` function.
  - The `PDB` module now defines the `query_alphafolddb` and `download_alphafold_structure`
    functions to query the *AlphaFold Protein Structure Database* and download the
    predicted structures.
  - This version solves a bug when reading MSA files with `|` in the sequence names.
  - MIToS is now using `Format.jl` instead of `Formatting.jl`.

### Changes from v2.14.1 to v2.15.0

  - The `MSA` module now exports the `rename_sequences!` and `rename_sequences` functions to
    rename the sequences of an MSA object.

### Changes from v2.14.0 to v2.14.1

  - The `modelled_sequences` function now returns only the selected chains, therefore avoid
    the inclusion of empty sequences in the output.

### Changes from v2.13.1 to v2.14.0

  - The `MSA` now defines `join` for MSA objects, allowing to join or merge two
    `AnnotationMultipleSequenceAlignment` objects based on a list of matching sequences or columns.

  - The `MSA` module now defines `hcat` and `vcat` for MSA objects, taking care of sequence
    and column names, and MSA annotations.
  - The `MSA` now exports the `sequencename_iterator` and `columnname_iterator` functions to
    return an iterator over the sequence or column names of an MSA.
  - The `MSA` now exports the `sequence_index` and `column_index` functions to return the
    integer position of a sequence or column name in an MSA.
  - `merge` and `merge!` are now defined for `Annotations` objects in the `MSA` module.

### Changes from v2.13.0 to v2.13.1

  - The `PDB` module can now parse the 66-character width columns of the PDB files created
    by *Foldseek*. These structures contain only the alpha carbons and do not have the column
    determining the element symbol.

### Changes from v2.12.0 to v2.13.0

  - The `PDB` module now includes the `modelled_sequences` function, allowing extraction of
    protein sequences from a specified structure.

  - The `PDB` module exports the `is_aminoacid` function to determine whether
    a `PDBResidue` represents an amino acid residue. This function is utilized by
    the `modelled_sequences` function.
  - The `Utils` module now exports the `THREE2ONE` constant, which is a dictionary mapping
    three-letter amino acid residue codes to their corresponding one-letter codes.

### Changes from v2.11.1 to v2.12.0

  - The `downloadsifts` function now downloads the SIFTS files from the PDBe HTTPS server
    instead of the previous FTP server. This improves error handling during the download
    process, making it more robust by relying on the `download_file` function. If you prefer
    the previous behavior, you can set the new keyword argument `source` to `"ftp"`.

  - It resolves an issue with the representation of Multiple Sequence Alignments and
    ContingencyTables in the `show` methods by always using explicit MIME types.
  - *[Breaking change]* The `show` methods that accept only two elements without an explicit
    MIME type are now deprecated.

### Changes from v2.11.0 to v2.11.1

  - MIToS now checks the magic number of gzip files immediately after download. If
    the gzip file does not have the correct header, MIToS will attempt to download
    it again. In Julia versions below 1.2, it will retry the download once. In
    Julia 1.2 or higher, it will retry the download five times, using an
    ExponentialBackOff.

### Changes from v2.10.0 to v2.11.0

  - *[breaking change]* `getCA` returns `missing` if a `PDBResidue` has no CA atom
    (before it was an `AssertionError`).

### Changes from v2.9.0 to v2.10.0

  - *[breaking change]* `downloadsifts` now uses `Base.download` instead of `download_file` as HTTP (1.7 or lower) doesn't support FTP. Because of that, it doesn't accept keywords argument as `download_file` besides `filename`.

  - MIToS now supports HTTP 1.0 and has migrated from using `HTTP.request` to using `HTTP.download` for `MIToS.Utils.download_file` dropping support on HTTP 0.8. Thanks, @kool7d!
  - The `downloadpfam` function now uses the InterPro API, as the [Pfam website has been discontinued](https://xfam.wordpress.com/2022/08/04/pfam-website-decommission/). Thanks, @timholy!
  - The `downloadpfam` function now has an `alignment` keyword argument for choosing which Pfam alignment download. The options are `"full"` (the default), `"seed"` and `"uniprot"`.
  - MIToS switched to GitHub Actions for CI. Thanks, @timholy!

### Changes from v2.8.6 to v2.9.0

  - New `matches` keyword argument in the `superimpose` function to determine the residues to be aligned. Thanks, @timholy!

### Changes from v2.8.1 to v2.8.6

  - You can pass keyword arguments from `downloadsifts` to `download_file`.

### Changes from v2.8.1 to v2.8.5

  - Fix bugs when concatenating concatenated MSAs using `hcat`.

### Changes from v2.8.1 to v2.8.4

  - Ensure that `gaussdca` use the correct project file.

### Changes from v2.8.1 to v2.8.3

  - Increase `PairwiseListMatrices` required version.

  - Fix bugs when concatenating concatenated MSAs using `hcat`.

### Changes from v2.8.0 to v2.8.1

Fix bug when `read`ing `hcat` generated MSA in `Stockholm` format.

### Changes from v2.7.0 to v2.8.0

Multiple bug fixes and improvements related to `getindex` and `hcat`.

  - *[breaking change]* MSA `getindex` can now change the order of the columns
    in an `AnnotatedMultipleSequenceAlignment`.

  - *[breaking change]* `convert` to MSA and sequence objects is now deprecated;
    use the corresponding constructor.
  - `gethcatmapping` to get the mapping to the concatenated MSAs.

### Changes from v2.6.1 to v2.7.0

  - *[breaking change]* MSA `getindex` with `:` or arrays now return an object of
    the same type. The annotations of an `AnnotatedMultipleSequenceAlignment` are
    modified according to the selection.

  - *[breaking change]* MSA `getindex` can now change the order of the sequences
    in an `AnnotatedMultipleSequenceAlignment`.
  - It adds `hcat` support for MSA objects, taking care of the MSA annotations.

### Changes from v2.6.0 to v2.6.1

  - `download_file` and other `download...` functions now use the proxy settings
    declared with the `HTTP_PROXY` and `HTTP_PROXY` environment variables.

### Changes from v2.5.0 to v2.6.0

  - The RESTful API of PDB has changed, and the Legacy Fetch API Web Service was shut down on
    December 9th, 2020. To adapt to the new changes, `PDBMLHeader` has been deprecated, and the
    `downloadpdbheader` and `getpdbdescription` functions now return different objects.

### Changes from v2.4.0 to v2.5.0

MIToS v2.5.0 drops support for *Julia 0.7* and adds support for *Julia 1.5* and
includes several bug fixes.

  - `Cookbook` section added to the docs using [Literate](https://github.com/fredrikekre/Literate.jl)

  - The `SIFTS` module now includes the `dbSCOP2` and `dbSCOP2B` databases.
  - `siftsmapping` now returns an `OrderedDict` instead of a `Dict`.
  - `msacolumn2pdbresidue` now return an `OrderedDict` instead of a `Dict`.

### Changes from v2.3.0 to v2.4.0

MIToS v2.4 uses `Project.toml` and includes several bug fixes.

  - The `SIFTS` module includes the `dbEnsembl` database and `warn`s again about unused databases.

### Changes from v2.2.0 to v2.3.0

MIToS v2.3 requires Julia v0.7 or v1.0. This release drops Julia 0.6 support.

  - `Formatting.jl` is used in place of `Format.jl`.

  - `SIFTS.get` returns the desired object or `missing` instead of `Nullable`s.
  - `SIFTS` function doesn't `warn` about unused databases.

#### Julia 0.7/1.0 deprecations

  - `bits` was deprecated to `bitstring`.

  - `'` and `.'` are deprecated for alignments and sequences, use `transpose` or
    `permutedims` instead. `ctranspose` is not longer available for matrices of `Residue`s.

### Changes from v2.1.2 to v2.2

  - `PIR` `FileFormat` is included to read and write alignments in PIR/NBRF format.

  - `Utils.Format` was renamed to `Utils.FileFormat`.
  - `HTTP.jl` is used in place of `FTPClient.jl` and the deprecated `Requests.jl` in
    `Utils.download_file` to download files.
  - `Format.jl` is used in place of `Formatting.jl`.
  - Solve bug in the printing of matrices of `Residue`s using `FileFormat`s.

### Changes from v2.1.1 to v2.1.2

  - `FTPClient.jl` is used in `Utils.download_file` to download files from FTP.

  - `CodecZlib.jl` is used in place of `GZip.jl` speeding up the parsing of compressed files.
  - Improvements in MSA and PDB parsing speed.
  - Improvement in `MSA.percentidentity` speed.
  - `Information.gaussdca` now uses Julia's `serialize` and `deserialize` instead of `JLD`.
  - `ROCAnalysis.jl` is not longer a dependency and it's now used with `@require` from
    `Requires.jl`. To use the `AUC` function you need to do `using ROCAnalysis`.

### Changes from v2.1 to v2.1.1

  - The script `Conservation.jl` was added to measure residue conservation of MSA columns.

  - The script `SplitStockholm.jl` now has a progress bar thanks to Ellis Valentiner
    @ellisvalentiner.

### Changes from v2.0 to v2.1

MIToS v2.1 requires Julia v0.6. This release drops Julia 0.5 support.

  - `get_n_words(...` doesn't remove the last newline character, use `get_n_words(chomp(...`
    to get the previous behaviour.

### Changes from v1.2.3 to v2.0

**MIToS 2.0** is the first MIToS version with **Julia 0.5** support
(It drops Julia 0.4 support). The last Julia version introduces new awesome features like
native multi-threading support, fast anonymous functions, generator expressions and more.
Also, the Julia package ecosystem has grown. So, MIToS was slightly redesigned to take
advantage of the new Julia capabilities. As a consequence, this version introduces several
breaking changes and new features.

##### Utils module

  - `deleteitems!(vector::Vector, items)` is deprecated in favor of
    `filter!(x -> x ∉ items, vector)`.

  - `All` is used instead of MIToS 1.0 `"all"` or `"*"`, because it's possible to dispatch on it.

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

  - `Residue` is now encoded as `Int` instead of being encoded as `UInt8`, allowing faster
    indexation using `Int(res::Residue)`. More memory is used, since the residues are encoded
    using 32 or 64 bits instead of 8 bits.

  - `XAA` is now used to indicate unknown, ambiguous and non standard residues instead of `GAP`.
  - Conversions to and from `UInt8` aren't supported now.
  - More `Base` methods are extended to work with `Residue`: `bits`, `zero`, `one`
    and `isvalid`.
  - `empty(Annotations)` was deprecated, use `Annotations()` instead.
  - `msa["seq_name",:]` now returns a `NamedArray{Residue,1}` instead of an aligned sequence,
    use `getsequence(msa,"seqname")` to get an aligned sequence with annotations.
  - The `names` function was replaced by the `sequencenames` function. A `columnnames`
    function was also added.
  - Aligned sequences don't drop dimensions, so there are matrices instead of vectors. You can
    use `vec(...)` or `squeeze(...,1)` to get a vector instead of the matrix.
  - Indexing MSA objects with only one string is deprecated, use `msa["seqname",:]` instead
    of `msa["seqname"]`.
  - `empty!` doesn't take MSA objects anymore.
  - `asciisequence` was replaced by `stringsequence`.
  - `deletenotalphabetsequences` and the parse/read keyword argument `checkalphabet` are
    deprecated since MIToS 2.0 uses Residue('X') to represent residues outside the alphabet. You
    can use `filtersequences!(msa, vec(mapslices(seq -> !in(XAA, seq), msa, 2)))` to delete
    sequences with unknown, ambiguous or non standard residues.
  - `parse`/`read` and MSA file returns an `AnnotatedMultipleSequenceAlignment` by default.
  - `shuffle_...columnwise!` and `shuffle_...sequencewise!` functions were deprecated in
    favor of `shuffle!` and `shuffle` functions.
  - `SequenceClusters` was renamed to `Clusters`.
  - Residue alphabet types were added. All alphabet types are subtypes of `ResidueAlphabet`.
    In particular, three types are exported: `GappedAlphabet`, `UngappedAlphabet` and
    `ReducedAlphabet`. The last type allows the creation of custom reduced alphabets.
  - In order to keep the sequence name, `AlignedSequence` and `AnnotatedAlignedSequence` are
    now matrices instead of vectors.

##### PDB module

  - The keyword argument `format` of `downloadpdb` should be a type (`PDBFile` or `PDBML`)
    instead of a string (`pdb` or `xml`) as in MIToS 1.0.

  - `read` and `parse` now has the `occupancyfilter` keyword argument.
  - `read` and `parse` now has the `label` keyword argument for `PDBML` files.
  - `residues`, `àtoms` and similiar functions don't take vectors or sets anymore. Use an
    anonymous function instead, e.g.: `x -> x in set_of_residue_numbers`.
  - The functions `isresidue`, `isatom` and `residuepairsmatrix` were added.

##### SIFTS module

  - The `get` function has a more complex signature for `SIFTSResidue`s to make simpler
    the access of data.

  - `find`, `filter` and `filter` now takes a database type as a third parameter when a vector
    of `SIFTSResidue`s is the second parameter. It allows to use a function that directly
    operates over the database type if it's available.
  - `SIFTSResidue`s now also store secondary structure data in the `sscode` and `ssname` fields.

##### Information module

  - `ResidueProbability` and `ResidueCount` were deprecated in favor of `ContingencyTable`.
    `Probabilities` and `Counts` were added as wrappers of `ContingencyTable` to allow dispach
    in a some functions, e.g. `entropy`.

  - The last parameter of contingency tables is now a subtype of `ResidueAlphabet` instead
    of a `Bool`, i.e.: `UngappedAlphabet`, `GappedAlphabet` or `ReducedAlphabet`.
  - Creation of empty contingecy tables chaged.
    e.g. `zeros(ResidueProbability{Float64, 2, false})` changed to
    `ContingencyTable(Float64, Val{2}, UngappedAlphabet())` and
    `ResidueProbability{Float64, 2, false}()` changed to
    `ContingencyTable{Float64, 2, UngappedAlphabet}(UngappedAlphabet())`.
  - `count!` and `probabilities!` signatures changed. The first argument is alway a
    `ContingencyTable`, the second positional argument a clustering weight object
    (use `NoClustering()` to skip it), the third positional argument is a pseudocount object
    (use `NoPseudocount()` to avoid the use of pseudocounts) and `probabilities!` takes also a
    `Pseudofrequencies` object (use `NoPseudofrequencies()` to avoid pseudofrequencies). The
    last positional arguments are the vector of residues used to fill the contingency table.
  - `count` and `probabilities` now takes the sequences as only positional arguments. The
    output is always a table of `Float64`. Both functions take the keyword arguments
    `alphabet`, `weights` and `pseudocounts`. `probabilities` also has a `pseudofrequencies`
    keyword argument.
  - `apply_pseudofrequencies!` changed its signature. Now it takes a `ContingencyTable` and
    a `Pseudofrequencies` object.
  - The function `blosum_pseudofrequencies!` was deprecated in favor of introducing a
    `BLOSUM_Pseudofrequencies` type as subtype of `Pseudofrequencies` to be used in
    `probabilities`, `probabilities!` and `apply_pseudofrequencies!`.
  - Because higher-order function are fast in Julia 0.5, measure types
    (i.e. subtypes of `AbstractMeasure`) were deprecated in favor of functions. In particular,
    `MutualInformation` was replaced with the `mutual_information` function,
    `MutualInformationOverEntropy` was replaced with `normalized_mutual_information`,
    `KullbackLeibler` was replaced with `kullback_leibler` and `Entropy` was replaced with
    `entropy`.
  - The functions `estimate`, `estimate_on_marginal` , `estimateincolumns` and
    `estimateinsequences` were deprecated because measure types are not longer used.
  - `estimate_on_marginal(Entropy...` was deprecated in favor of the `marginal_entropy`
    function.
  - `estimateincolumns` and `estimateinsequences` were deprecated in favor of `mapcolfreq!`,
    `mapseqfreq!`, `mapcolpairfreq!` and `mapseqpairfreq`.
  - Keyword argument `usegaps` is deprecated in `buslje09` and `BLMI` in favor of `alphabet`.
  - `cumulative` function was added to calculate cumulative MI (cMI).

* * *

### Changes from v1.1 to v1.2.2

  - `using Plots` to use `plot` with `AbstractVector{PDBResidue}` to visualize coordinates
    of the C alpha of each residue.

  - Re-exports `swap!` from **IndexedArrays.jl**.
  - *[breaking change]* **Distances.jl** now uses `--inter` instead of `--intra`.
  - *docs* and *cookbook* are now in [MIToSDocumentation](https://github.com/diegozea/MIToSDocumentation)

* * *

### Changes from v1.0 to v1.1

  - **RecipesBase** is used to generate plot recipes for MIToS’ objects. MSA objects can be
    visualized `using Plots` (thanks to Thomas Breloff @tbreloff ).

  - Functions to perform structural superimposition were added to the `PDB` module
    (thanks to Jorge Fernández de Cossío Díaz @cosio ) : `center!`, `kabsch`, `rmsd`.
  - The `PDB` module adds the following functions to make easier structural comparison:
    `getCA`, `CAmatrix`, `coordinatesmatrix`, `centeredcoordinates`, `centeredresidues`,
    `change_coordinates`, `superimpose`, `mean_coordinates` and `rmsf`.
  - When PDB or PDBML files are being parsed, It’s possible to indicate if only atoms with
    the best occupancy should be loaded (`occupancyfilter=true`, `false` by default).
  - When `PDBML` files are being parsed, is possible to used the new `label` keyword argument
    to indicate if "auth" (`false`) or "label" (`true`) attributes should be used.
  - `bestoccupancy!` was deprecated in favor of `bestoccupancy`.
  - The `MSA` module export the function `percentsimilarity` to calculate the similarity
    percent between aligned sequences.
  - `msacolumn2pdbresidue` has two new keyword arguments, `strict` and `checkpdbname`, to
    perform extra tests during the mapping between PDB and MSA residues.
  - `msacolumn2pdbresidue` has a new `missings` keyword argument to indicate if missing
    residues should be included in the mapping (default: `true`).
  - The `MSA` now exports the `residue2three` and `three2residue` function to convert
    `Residue`s to and from their three letter names.
  - The `MSA` module now exports `sequencepairsmatrix`, `columnpairsmatrix`, `columnlabels`,
    and `sequencelabels` to help in the construction of matrices for MSA sequences or columns
    pairwise comparisons.
  - The `Information` module, if `GaussDCA` is installed, allows to call its `gDCA` function
    from MIToS through the `gaussdca` function.
  - The `Information` module now exports the `KullbackLeibler` measure.
  - Now is possible to `print` and `write` `PDBResidue`s as `PDBFile`s.
  - The function `proximitymean` now has a keyword argument `include` to indicate if the
    residue score should be included in the mean.
  - The module `Scripts` inside the `Utils` module has a new function `readorparse` to help
    parsing `STDIN` in MIToS’ scripts.

**MIToS v1.1** also includes several **bug fixes**, some **performance improvements** and a
more complete **documentation**.

* * *

### Changes from v0.1 to v1.0

  - `Pfam` module for working with *Pfam* alignments and useful parameter optimization
    functions (i.e. `AUC`).

  - *[breaking change]* The `Clustering` module was deleted and its functions moved to the
    `MSA` module.
  - `MSA` uses `ClusteringResult` from the `Clustering.jl` package instead of `AbstractClusters`.
    
      + `Clusters` was renamed to `SequenceClusters`
    
      + `MSA` adds the `counts` and `assignments` functions from the `Clustering.jl` interface.
      + *[breaking change]* The `getnclusters` function is now `nclusters` in the `Clutering` module.
  - *[breaking change]* All the MSA `...percentage` functions were renamed to `...fraction`
    and `percent...` functions now return real percentages (not fractions) values.
    Functions taking identity thresholds, now also take real percentages
    (values between 0.0 and 100.0).
  - *[breaking change]* Script command line arguments changed to: define the number of
    workers, use STDIN and STDOUT (pipelines), get better output names, use real flag arguments.
  - `InformationMeasure` renamed to `AbstractMeasure`.
  - New functions added to `MSA` module.
    
      + `annotations`, `names`.
    
      + `meanpercentidentity` allows fast estimation of the mean percent identity between the sequences of a MSA.
  - New function and type added to `Information` module.
    
      + `cumulative` to calculate cMI (cumulative mutual information) and similar cumulative scores.
    
      + `KullbackLeibler` to estimate conservation.
  - `proximitymean` is defined in the `PDB` module to calculate pMI
    (proximity mutual information) and other proximity scores.
  - `contact` and `distance` have a vectorized form to create contact/distance maps.
  - `NCol` file annotation with the number of columns in the original MSA.
  - `BLMI` has `lambda` as a keyword argument for using additive smoothing.
  - `BLMI` and `buslje09` accepts `samples=0` to avoid the Z score estimation.
  - `read`/`parse` added the keyword argument `checkalphabet` for deleting sequences with non
    standard amino acids.
  - `read`/`parse` added the keyword argument `keepinserts` for keep insert columns
    (It creates an `Aligned` column annotation).

**MIToS v1.0** also includes several **bug fixes** and a more complete **documentation**.
