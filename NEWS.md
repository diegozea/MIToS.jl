## MIToS.jl Release Notes

### Changes from v0.1 to v0.2

* *[breaking change]* Percent identity functions were moved from `Clutering` to `MSA` module.

* *[breaking change]* The `getnclusters` function is now `nclusters` in the `Clutering` module.

* *[breaking change]* All the MSA `...percentage` functions were renamed to `...fraction`.

* New functions added to `MSA` module

  * `annotations`, `names`

  * `meanpercentidentity` allows fast estimation of the mean percent identity between the sequences of a MSA.

* `NCol` file annotation with the number of columns in the original MSA.

* `Pfam` module for working with *Pfam* alignments and useful parameter optimization functions.

* `BLMI` has `lambda` as a keyword argument for using additive smoothing.

* `read`/`parse` added the keyword argument `checkalphabet` for deleting sequences with non standard amino acids.

* `BLMI` and `buslje09` accepts `samples=0` to avoid the Z score estimation.

**MIToS v0.2** also includes several **bug fixes** and a more complete **documentation**.
