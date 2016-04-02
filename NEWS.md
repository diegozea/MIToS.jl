## MIToS.jl Release Notes

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
