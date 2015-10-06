# MIToS
### Mutual Information Tools for protein Sequence analysis

Linux, OSX: [![Build Status](https://travis-ci.org/diegozea/MIToS.jl.svg?branch=master)](https://travis-ci.org/diegozea/MIToS.jl)

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/h6o72b5dtdeto336/branch/master?svg=true)](https://ci.appveyor.com/project/diegozea/mitos-jl/branch/master)

Code Coverage: [![Coverage Status](https://coveralls.io/repos/diegozea/MIToS.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/diegozea/MIToS.jl?branch=master) [![codecov.io](http://codecov.io/github/diegozea/MIToS.jl/coverage.svg?branch=master)](http://codecov.io/github/diegozea/MIToS.jl?branch=master)

This implementation of Mutual Information (MI) analysis and several useful tools for Multiple Sequence Alignments (MSA)s and PDB structures management in Julia language. MI is useful for determining covariation between positions in a MSA. Mutual information derived scores can be useful for determining structural contacts and functional sites on proteins [1,2].

MIToS starting point was an improvement of the algorithm published by Buslje et. al. [1]. A BLOSUM62-based pseudocount strategy, similar to Altschul et. al. [3], was implemented for getting a better performance with few number of sequences in the MSA. **MIToS** implements all the necessary tools for using, developing and testing scores based on MI in different modules:   

* **MSA** defines multiple functions and types for dealing with MSAs:  
  * `AnnotatedMultipleSequenceAlignment` is a type for saves a MSA and its `Annotations`, useful for working with MSAs on `Stockholm`/Pfam format. 
  * `Annotations` can store the sequence and column mapping after operations like `gapstrip!` or `adjustreference!`.
  * Read and write `FASTA`, `Stockholm` or `Raw` formats.
  * Functions for shuffling the MSA: `shuffle_...`
* **Clustering** defines a `Clusters` type for clustering sequences from a MSA and implements Hobohm I [4] for sequence weighting: `hobohmI`  
* **PDB** defines fuctions for parsing and working with `PDBFile` and `PDBML` formats:  
  * Defines the types: `PDBResidue`, `PDBResidueIdentifier`, `PDBAtom`, `Coordinates`
  * Functions for estimation of `distance` and `contact` between `PDBResidue`s, also another type of interactions like: `vanderwaals`, `ionic`, etc.
  * Macros and functions for getting residues or atoms from a list/vector of residues, i.e.: `@residues`
* **SIFTS** has functions for download and parse per PDB SIFTS XML files:  
  * `DataBase` and subtypes `db...` (i.e. `dbUniProt`) for a residue level mapping between databases.
  * `siftsmapping` function for getting a easy-to-use `Dict` mapping. 
  * `SIFTSResidue` type for a low level interface to the SIFTS mapping.
* **Information** has functions and types for measuring information content:  
  * BLOSUM62 probabilities: `BLOSUM62_Pi` and `BLOSUM62_Pij`
  * `ResidueContingencyTables` for counting or store the probabilities of `Residues` (using the 20 residues BLOSUM62 alphabet and allows to count gaps)
  * `InformationMeasure`s:
    * `Entropy`
    * `MutualInformation`
    * `MutualInformationOverEntropy`
  * Corrections to the MI values for co-evolution estimation:
    * `APC!` for MIp
    * `buslje09` for the Z score and MIp from Buslje et. al. [1]
* **Utils** for common utils functions and types for MIToS.

#### Scripts   

**MIToS** implements several useful scripts for command line execution (without requiring Julia coding):  
  
* **Buslje09.jl** : Calculates a Z score and a corrected MI/MIp as described on Buslje et. al. [1].
* **DownloadPDB.jl** : Download gzipped files from PDB.
* **Distances.jl** : Calculates residues distance from a PDB file.
* **SplitStockholm.jl** : Splits a Stockholm file with multiple alignments into one compressed file per MSA
* **AlignedColumns.jl** : Creates a Stockholm file with the aligned columns from a Pfam Stockholm file (insertions are deleted) saving the mapping for the sequences (residue number in UniProt) and the columns in the original MSA.
* **MSADescription.jl** : Calclulates from a Stockholm file the number of columns, sequences and clusters after Hobohm I clustering at 62% identity. Also gives the mean, standard deviation and quantiles of: sequence coverage of the MSA, gap percentage.

#### References  

1. Buslje, C. M., Santos, J., Delfino, J. M., & Nielsen, M. **Correction for phylogeny, small number of observations and data redundancy improves the identification of coevolving amino acid pairs using mutual information.** *Bioinformatics 2009*, 25(9), 1125-1131.  
2. Buslje, C. M., Teppa, E., Di Doménico, T., Delfino, J. M., & Nielsen, M. **Networks of high mutual information define the structural proximity of catalytic sites: implications for catalytic residue identification.** *PLoS Comput Biol 2010*, 6(11), e1000978-e1000978.
3. Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J., Zhang, Z., Miller, W., & Lipman, D. J. **Gapped BLAST and PSI-BLAST: a new generation of protein database search programs.** *Nucleic acids research 1997*, 25(17), 3389-3402.
4. Hobohm, U., Scharf, M., Schneider, R., & Sander, C. **Selection of representative protein data sets.** *Protein Science 1992*, 1(3), 409-417.


#### Structural Bioinformatics Unit  
[![FIL](http://mistic.leloir.org.ar/imgs/logo_horizontal.png)](http://www.leloir.org.ar/)
