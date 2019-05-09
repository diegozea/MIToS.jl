[![MIToS](http://diegozea.github.io/MIToS.jl/latest/assets/MIToS_logo.png)](http://diegozea.github.io/MIToS.jl/)
### Mutual Information Tools for protein Sequence analysis

Julia 0.7 & 1.0: **MIToS 2.4.0**  

Linux, OSX:
[![Build Status](https://travis-ci.org/diegozea/MIToS.jl.svg?branch=master)](https://travis-ci.org/diegozea/MIToS.jl)  
Windows:
[![Build status](https://ci.appveyor.com/api/projects/status/h6o72b5dtdeto336/branch/master?svg=true)](https://ci.appveyor.com/project/diegozea/mitos-jl/branch/master)  

Code Coverage:
[![Coverage Status](https://coveralls.io/repos/diegozea/MIToS.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/diegozea/MIToS.jl?branch=master) [![codecov.io](http://codecov.io/github/diegozea/MIToS.jl/coverage.svg?branch=master)](http://codecov.io/github/diegozea/MIToS.jl?branch=master)

**NOTE:** [Julia has reached a stable 1.0 version](https://julialang.org/blog/2018/08/one-point-zero) and **MIToS v2.4** supports it. You need to do `using Pkg; Pkg.add("MIToS")` to **install** MIToS on Julia 1.0 or `using Pkg; Pkg.update()` to update your installed version.  

Some breaking changes were introduced in MIToS v2.3. See the [NEWS.md](https://github.com/diegozea/MIToS.jl/blob/master/NEWS.md)
file and the [documentation](https://diegozea.github.io/MIToS.jl/stable) to migrate code from an old version
of MIToS. If you need more help to migrate code from MIToS 1.0 in Julia 0.4 to MIToS 2.4 in Julia 1.0, you can
write a mail to diegozea at gmail dot com asking for assistance.  

**DOCUMENTATION**:  
Documentation for [MIToS 1.0 in Julia 0.4](https://diegozea.github.io/mitosghpage-legacy/)  
Documentation for MIToS 2.0 or greater in Julia 0.5 or greater: [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://diegozea.github.io/MIToS.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://diegozea.github.io/MIToS.jl/latest)  

MIToS is an environment for Mutual Information (MI) analysis and implements several useful
tools for Multiple Sequence Alignments (MSAs) and PDB structures management in the Julia
language [1]. MI allows determining covariation between positions in a MSA. MI derived scores
are good predictors of residue contacts and functional sites in proteins [2,3].

MIToS starting point was an improvement of the algorithm published by Buslje et. al. [2]. A BLOSUM62-based pseudocount strategy, similar to Altschul et. al.[4], was implemented for
a better performance in the range of MSAs with low number of sequences. **MIToS** offers
all the necessary tools for using, developing and testing MI based scores, in different
modules.

#### Modules
MIToS tools are separated on different modules, related to different tasks.
- **MSA** This module defines multiple functions and types for dealing with MSAs and
their annotations. It also includes facilities for sequence clustering.
- **PDB** This module defines types and methods to work with protein structures from PDB.
- **SIFTS** This module allows access to SIFTS residue-level mapping of UniProt, Pfam and
other databases with PDB entries.
- **Information** This module defines residue contingency tables and methods on them
to estimate information measure from MSAs. It includes functions to estimate corrected
mutual information (ZMIp, ZBLMIp) between MSA columns.
- **Pfam**
This module use the previous modules to work with Pfam MSAs. It also has useful parameter
optimization functions to be used with Pfam alignments.
- **Utils** MIToS has also an Utils module with common utils functions and types used
in this package.

#### Scripts

**MIToS** implements several useful scripts for command line execution
(without requiring Julia coding):

* **Buslje09.jl** : Calculates the corrected MI/MIp described on Buslje et. al. 2009 [2].
* **BLMI.jl** : Calculates corrected mutual information using BLOSUM62 based-pseudocounts.
* **Conservation.jl** : Calculates the Shannon entropy and the Kullback-Leibler divergence
of each MSA column.
* **DownloadPDB.jl** : Downloads gzipped files from PDB.
* **Distances.jl** : Calculates residues distances in a PDB file.
* **SplitStockholm.jl** : Splits a Stockholm file with multiple alignments into one
compressed file per MSA
* **AlignedColumns.jl** : Creates a Stockholm file with the aligned columns from a Pfam
Stockholm file (insertions are deleted) saving the mapping (residue number in UniProt)
and the columns in the original MSA.
* **PercentIdentity.jl** : Calculates the percentage identity between all the sequences
of an MSA and saves mean, median, minimum, etc.
* **MSADescription.jl** : Calculates the number of columns, sequences and clusters after
Hobohm I clustering at 62% identity given a stockholm file as input [5]. It also gives the
percent indentity mean and mean, standard deviation and quantiles of: sequence coverage of
the MSA and gap percentage.

#### Citation  

If you use MIToS, please cite:

Diego J. Zea, Diego Anfossi, Morten Nielsen, Cristina Marino-Buslje; **MIToS.jl: mutual information tools for protein sequence analysis in the Julia language**, Bioinformatics, Volume 33, Issue 4, 15 February 2017, Pages 564–565, [https://doi.org/10.1093/bioinformatics/btw646](https://doi.org/10.1093/bioinformatics/btw646)

#### References

1. Zea, Diego Javier, et al. "MIToS. jl: mutual information tools for protein sequence
analysis in the Julia language." Bioinformatics 33, no. 4 (2016): 564-565.
2. Buslje, Cristina Marino, et al. "Correction for phylogeny, small number of
observations and data redundancy improves the identification of coevolving amino acid
pairs using mutual information." Bioinformatics 25.9 (2009): 1125-1131.
3. Buslje, Cristina Marino, et al. "Networks of high mutual information define the
structural proximity of catalytic sites: implications for catalytic residue
identification." PLoS Comput Biol 6.11 (2010): e1000978.
4. Altschul, Stephen F., et al. "Gapped BLAST and PSI-BLAST: a new generation of protein
database search programs." Nucleic acids research 25.17 (1997): 3389-3402.
5. Hobohm, Uwe, et al. "Selection of representative protein data sets." Protein Science
1.3 (1992): 409-417.

#### Structural Bioinformatics Unit
[![FIL](http://mistic.leloir.org.ar/imgs/logo_horizontal.png)](http://www.leloir.org.ar/)
