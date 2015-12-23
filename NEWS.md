## MIToS.jl Release Notes

### Changes from v0.1 to v0.2

* `BLMI` has `lambda` as a keyword argument for using additive smoothing.

* `read`/`parse` added the keyword argument `checkalphabet` for deleting sequences with non standard amino acids.

* `BLMI` and `buslje09` accepts `samples=0` to avoid the Z score estimation.

* `meanpercentidentity` allows fast estimation of the mean percent identity between the sequences of a MSA.

**MIToS v0.2** also includes several bug fixes.