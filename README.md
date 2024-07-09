[![MIToS](http://diegozea.github.io/MIToS.jl/latest/assets/MIToS_logo.png)](http://diegozea.github.io/MIToS.jl/)
## 🐉 MIToS: Mutual Information Tools for protein Sequence analysis

*A Julia Package to Analyze Protein Sequences, Structures, and Evolutionary Information*

<br>

**DOCUMENTATION:** [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://diegozea.github.io/MIToS.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://diegozea.github.io/MIToS.jl/latest)  

Linux, OSX & Windows: [![Status](https://github.com/diegozea/MIToS.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/diegozea/MIToS.jl/actions?query=workflow%3A%22CI%22+branch%3Amaster) Code Coverage:
[![Coverage Status](https://coveralls.io/repos/diegozea/MIToS.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/diegozea/MIToS.jl?branch=master) [![codecov.io](http://codecov.io/github/diegozea/MIToS.jl/coverage.svg?branch=master)](http://codecov.io/github/diegozea/MIToS.jl?branch=master)

> **NOTE:**  Some **breaking changes** were introduced between **MIToS 2.15** and **MIToS 3.0**, inclusive. See the [NEWS.md](https://github.com/diegozea/MIToS.jl/blob/master/NEWS.md) file to migrate code from an old version
of MIToS. Most breaking changes will show a deprecation warning with a hint on how to perform the migration. If you need more help migrating code towards MIToS v3, you can write an email to diegozea at gmail dot com asking for assistance.

MIToS provides a comprehensive suite of tools for the analysis of protein sequences and structures.
It allows working with **Multiple Sequence Alignments (MSAs)** to obtain evolutionary information in the Julia language [1].
In particular, it eases the analysis of coevoling position in an MSA using **Mutual Information (MI)**, a measure of covariation.
MI-derived scores are good predictors of inter-residue contacts in a protein structure and functional sites in proteins [2,3].
To allow such analysis, MIToS also implements several useful tools for working with protein structures, such as those available in the **Protein Data Bank (PDB)** or predicted by AlphaFold 2.

MIToS starting point was an improvement of the algorithm published by Buslje et al. [2]. 
A BLOSUM62-based pseudo-count strategy, similar to Altschul et al.[4], , was implemented to improve performance in the range of MSAs with a low number of sequences [1]. 
**MIToS** offers all the tools for using, developing, and testing MI-based scores—in fact, any measure based on reside frequencies in an MSA—in different modules.

### Modules
MIToS tools are separated into different modules for different tasks.
- **MSA** This module defines multiple functions and types for dealing with MSAs and
their annotations. It also includes facilities for sequence clustering and shuffling, among others.
- **PDB** This module defines types and methods to work with protein structures from
different sources, such as PDB or AlphaFold DB. It includes functions to superpose structures,
measure the distance between residues, and much more.
- **Information** This module defines residue contingency tables and methods 
to estimate information measures from MSAs. It includes functions to estimate corrected
mutual information (ZMIp, ZBLMIp) between MSA columns, as well as conservation estimations
using Shannon entropy and the Kullback-Leibler divergence.
- **SIFTS** This module allows access to SIFTS residue-level mapping of UniProt, Pfam, and
other databases with PDB entries.
- **Pfam** This module uses the previous modules to work with Pfam MSAs. It also offers
useful functions for parameter optimization using Pfam alignments.
- **Utils** It exports common utils functions and types used in different modules of this package.

### Installation

To install MIToS, you need to execute the following code in Julia:

```julia
using Pkg; Pkg.add("MIToS")
```
 
To update your installed version, you can execute:

```julia
using Pkg; Pkg.update("MIToS")`
```

### Scripts
The [MIToS_Scripts](https://github.com/MIToSOrg/MIToS_Scripts.jl) package offers a set of easy-to-use scripts to access some functionalities MIToS offers from the terminal. These scripts are designed for researchers familiar with command-line interfaces (CLI) but without experience coding in Julia. The available scripts include:

* **Buslje09.jl**: Calculates corrected Mutual Information (MI/MIp) based on Buslje et al., 2009.
* **BLMI.jl**: Computes corrected mutual information using BLOSUM62-based pseudo-counts, as described in the MIToS publication [1].
* **Conservation.jl**: Calculates Shannon entropy and Kullback-Leibler divergence for each MSA column.
* **Distances.jl**: Computes inter-residue distances in a PDB file.
* **PercentIdentity.jl**: Calculates the percentage identity between all sequences in an MSA and provides statistical summaries.
* **MSADescription.jl**: Provides statistics for a given Stockholm file, including clustering information and sequence coverage.

This list is not exhaustive; more scripts are available in the [MIToS_Scripts.jl repository](https://github.com/MIToSOrg/MIToS_Scripts.jl). Visit the repository for more details and to access these scripts.

### Order versions
MIToS 3.0 requires Julia 1.9 or higher. It is recommended that you use these versions to get the best experience coding with Julia and MIToS.
If you need to use MIToS in a Julia version lower than 1.0, you will need to look at the [older MIToS v1 documentation](https://diegozea.github.io/mitosghpage-legacy/).

### Citation  
If you use MIToS, please cite:

Diego J. Zea, Diego Anfossi, Morten Nielsen, Cristina Marino-Buslje; **MIToS.jl: mutual information tools for protein sequence analysis in the Julia language**, Bioinformatics, Volume 33, Issue 4, 15 February 2017, Pages 564–565, [https://doi.org/10.1093/bioinformatics/btw646](https://doi.org/10.1093/bioinformatics/btw646)

### References

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

### Structural Bioinformatics Unit
[![FIL](http://mistic.leloir.org.ar/imgs/logo_horizontal.png)](http://www.leloir.org.ar/)
