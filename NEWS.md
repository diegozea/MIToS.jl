## MIToS.jl Release Notes

### Changes from v0.1 to v0.2

* *[breaking change]* The `Clustering` module was deleted and its functions moved to the `MSA` module.

  * `MSA` uses `ClusteringResult` from the `Clustering.jl` package instead of `AbstractClusters`.

  * `Clusters` was renamed to `SequenceClusters`

  * `MSA` adds the `counts` and `assignments` functions from the `Clustering.jl` interface.

* *[breaking change]* The `getnclusters` function is now `nclusters` in the `Clutering` module.

* *[breaking change]* All the MSA `...percentage` functions were renamed to `...fraction` and `percent...` functions now return real percentages (not fractions) values.
Functions taking identity thresholds, now also take real percentages (values between 0.0 and 100.0).

* `InformationMeasure` renamed to `AbstractMeasure`

* New functions added to `MSA` module

  * `annotations`, `names`

  * `meanpercentidentity` allows fast estimation of the mean percent identity between the sequences of a MSA.

* `NCol` file annotation with the number of columns in the original MSA.

* `Pfam` module for working with *Pfam* alignments and useful parameter optimization functions.

* `BLMI` has `lambda` as a keyword argument for using additive smoothing.

* `read`/`parse` added the keyword argument `checkalphabet` for deleting sequences with non standard amino acids.

* `BLMI` and `buslje09` accepts `samples=0` to avoid the Z score estimation.

**MIToS v0.2** also includes several **bug fixes** and a more complete **documentation**.
