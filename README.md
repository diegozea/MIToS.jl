[![MIToS](http://diegozea.github.io/MIToS.jl/assets/MIToS_logo.png)](http://diegozea.github.io/MIToS.jl/)
### Mutual Information Tools for protein Sequence analysis

Julia 0.4: [![MIToS](http://pkg.julialang.org/badges/MIToS_0.4.svg)](http://pkg.julialang.org/?pkg=MIToS)

Linux, OSX: [![Build Status](https://travis-ci.org/diegozea/MIToS.jl.svg?branch=master)](https://travis-ci.org/diegozea/MIToS.jl)

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/h6o72b5dtdeto336/branch/master?svg=true)](https://ci.appveyor.com/project/diegozea/mitos-jl/branch/master)

Code Coverage: [![Coverage Status](https://coveralls.io/repos/diegozea/MIToS.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/diegozea/MIToS.jl?branch=master) [![codecov.io](http://codecov.io/github/diegozea/MIToS.jl/coverage.svg?branch=master)](http://codecov.io/github/diegozea/MIToS.jl?branch=master)

[**Documentation**](http://mitos.leloir.org.ar) : [http://mitos.leloir.org.ar](http://mitos.leloir.org.ar)

MIToS is an environment for Mutual Information (MI) analysis and implements several useful tools for Multiple Sequence Alignments (MSAs) and PDB structures management in the Julia language. MI allows determining covariation between positions in a MSA. MI derived scores are good predictors of residue contacts and functional sites in proteins [1,2].

MIToS starting point was an improvement of the algorithm published by Buslje et. al. [1]. A BLOSUM62-based pseudocount strategy, similar to Altschul et. al. [3], was implemented for a better performance in the range of MSAs with low number of sequences. **MIToS** offers all the necessary tools for using, developing and testing MI based scores, in different modules.

#### Modules
MIToS tools are separated on different modules, related to different tasks.
- **MSA** This module defines multiple functions and types for dealing with MSAs and their annotations. It also includes facilities for sequence clustering.
- **PDB** This module defines types and methods to work with protein structures from PDB.
- **SIFTS** This module allows access to SIFTS residue-level mapping of UniProt, Pfam and other databases with PDB entries.
- **Information** This module defines residue contingency tables and methods on them to estimate information measure from MSAs. It includes functions to estimate corrected mutual information (ZMIp, ZBLMIp) between MSA columns.
- **Pfam** 
This module use the previous modules to work with Pfam MSAs. It also has useful parameter optimization functions to be used with Pfam alignments.
- **Utils** MIToS has also an Utils module with common utils functions and types used in this package.

#### Scripts

**MIToS** implements several useful scripts for command line execution (without requiring Julia coding):

* **Buslje09.jl** : Calculates the corrected MI/MIp described on Buslje et. al. 2009 [1].
* **BLMI.jl** : Calculates corrected mutual information using BLOSUM62 based-pseudocounts.
* **DownloadPDB.jl** : Downloads gzipped files from PDB.
* **Distances.jl** : Calculates residues distances in a PDB file.
* **SplitStockholm.jl** : Splits a Stockholm file with multiple alignments into one compressed file per MSA
* **AlignedColumns.jl** : Creates a Stockholm file with the aligned columns from a Pfam Stockholm file (insertions are deleted) saving the mapping (residue number in UniProt) and the columns in the original MSA.
* **PercentIdentity.jl** : Calculates the percentage identity between all the sequences of an MSA and saves mean, median, minimum, etc.
* **MSADescription.jl** : Calculates the number of columns, sequences and clusters after Hobohm I clustering at 62% identity given a stockholm file as imput. It also gives the percent indentity mean and mean, standard deviation and quantiles of: sequence coverage of the MSA and gap percentage.

#### References

1. Buslje, Cristina Marino, et al. "Correction for phylogeny, small number of observations and data redundancy improves the identification of coevolving amino acid pairs using mutual information." Bioinformatics 25.9 (2009): 1125-1131.
2. Buslje, Cristina Marino, et al. "Networks of high mutual information define the structural proximity of catalytic sites: implications for catalytic residue identification." PLoS Comput Biol 6.11 (2010): e1000978.
3. Altschul, Stephen F., et al. "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic acids research 25.17 (1997): 3389-3402.
4. Hobohm, Uwe, et al. "Selection of representative protein data sets." Protein Science 1.3 (1992): 409-417.
5. Dunn, Stanley D., Lindi M. Wahl, and Gregory B. Gloor. "Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction." Bioinformatics 24.3 (2008): 333-340.

#### Structural Bioinformatics Unit
[![FIL](http://mistic.leloir.org.ar/imgs/logo_horizontal.png)](http://www.leloir.org.ar/)
