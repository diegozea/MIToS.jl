## MIToS.jl Release Notes

### Changes from v1.1 to v1.2

* `using Plots` to use `plot` with `AbstractVector{PDBResidue}` to visualize coordinates of the C alpha of each residue

* *[breaking change]* **Distances.jl** now uses `--inter` instead of `--intra`

### Changes from v1.0 to v1.1

* **RecipesBase** is used to generate plot recipes for MIToS’ objects. MSA objects can be visualized `using Plots` (thanks to Thomas Breloff @tbreloff ).

* Functions to perform structural superimposition were added to the `PDB` module (thanks to Jorge Fernández de Cossío Díaz @cosio ) : `center!`, `kabsch`, `rmsd`.

* The `PDB` module adds the following functions to make easier structural comparison: `getCA`, `CAmatrix`, `coordinatesmatrix`, `centeredcoordinates`, `centeredresidues`, `change_coordinates`, `superimpose`, `mean_coordinates` and `rmsf`.

* When PDB or PDBML files are being parsed, It’s possible to indicate if only atoms with the best occupancy should be loaded (`occupancyfilter=true`, `false` by default).

* When `PDBML` files are being parsed, is possible to used the new `label` keyword argument to indicate if "auth" (`false`) or "label" (`true`) attributes should be used.

* `bestoccupancy!` was deprecated in favor of `bestoccupancy`.

* The `MSA` module export the function `percentsimilarity` to calculate the similarity percent between aligned sequences.

* `msacolumn2pdbresidue` has two new keyword arguments, `strict` and `checkpdbname`, to perform extra tests during the mapping between PDB and MSA residues.

* `msacolumn2pdbresidue` has a new `missings` keyword argument to indicate if missing residues should be included in the mapping (default: `true`).

* The `MSA` now exports the `residue2three` and `three2residue` function to convert `Residue`s to and from their three letter names.

* The `MSA` module now exports `sequencepairsmatrix`, `columnpairsmatrix`, `columnlabels`, and `sequencelabels` to help in the construction of matrices for MSA sequences or columns pairwise comparisons.

* The `Information` module, if `GaussDCA` is installed, allows to call its `gDCA` function from MIToS through the `gaussdca` function.

* The `Information` module now exports the `KullbackLeibler` measure.

* Now is possible to `print` and `write` `PDBResidue`s as `PDBFile`s.

* The function `proximitymean` now has a keyword argument `include` to indicate if the residue score should be included in the mean.

* The module `Scripts` inside the `Utils` module has a new function `readorparse` to help parsing `STDIN` in MIToS’ scripts.

**MIToS v1.1** also includes several **bug fixes**, some **performance improvements** and a more complete **documentation**.

### Changes from v0.1 to v1.0

* `Pfam` module for working with *Pfam* alignments and useful parameter optimization functions (i.e. `AUC`).

* *[breaking change]* The `Clustering` module was deleted and its functions moved to the `MSA` module.

  * `MSA` uses `ClusteringResult` from the `Clustering.jl` package instead of `AbstractClusters`.

  * `Clusters` was renamed to `SequenceClusters`

  * `MSA` adds the `counts` and `assignments` functions from the `Clustering.jl` interface.

* *[breaking change]* The `getnclusters` function is now `nclusters` in the `Clutering` module.

* *[breaking change]* All the MSA `...percentage` functions were renamed to `...fraction` and `percent...` functions now return real percentages (not fractions) values.
Functions taking identity thresholds, now also take real percentages (values between 0.0 and 100.0).

* *[breaking change]* Script command line arguments changed to: define the number of workers, use STDIN and STDOUT (pipelines), get better output names, use real flag arguments.

* `InformationMeasure` renamed to `AbstractMeasure`.

* New functions added to `MSA` module.

  * `annotations`, `names`.

  * `meanpercentidentity` allows fast estimation of the mean percent identity between the sequences of a MSA.

* New function and type added to `Information` module.

  * `cumulative` to calculate cMI (cumulative mutual information) and similar cumulative scores.

  * `KullbackLeibler` to estimate conservation.

* `proximitymean` is defined in the `PDB` module to calculate pMI (proximity mutual information) and other proximity scores.

* `contact` and `distance` have a vectorized form to create contact/distance maps.

* `NCol` file annotation with the number of columns in the original MSA.

* `BLMI` has `lambda` as a keyword argument for using additive smoothing.

* `BLMI` and `buslje09` accepts `samples=0` to avoid the Z score estimation.

* `read`/`parse` added the keyword argument `checkalphabet` for deleting sequences with non standard amino acids.

* `read`/`parse` added the keyword argument `keepinserts` for keep insert columns (It creates an `Aligned` column annotation).

**MIToS v1.0** also includes several **bug fixes** and a more complete **documentation**.
