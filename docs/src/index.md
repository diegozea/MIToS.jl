![MIToS.jl](./assets/MIToS_logo.png)  

*A Julia Package to Analyze Protein Sequences, Structures, and Evolutionary Information*

## Modules

MIToS tools are separated into different modules for different tasks.

- [MSA](@ref Module-MSA): This module defines multiple functions and types for dealing with 
Multiple Sequence Alignments (MSAs) and their annotations. It also includes facilities 
for sequence clustering and shuffling, among others.

- [PDB](@ref Module-PDB): This module defines types and methods to work with protein 
structures from different sources, such as the Protein Data Bank (PDB) or AlphaFold DB. 
It includes functions to superpose structures, measure the distance between residues, and much more.

- [Information](@ref Module-Information): This module defines residue contingency tables and 
methods on them to estimate information measures. This allow to measure evolutionary
information on MSAs positions. It includes functions to estimate corrected mutual 
information (ZMIp, ZBLMIp) between MSA columns, as well as conservation estimations using 
Shannon entropy and the Kullback-Leibler divergence.

- [SIFTS](@ref Module-SIFTS): This module allows access to SIFTS residue-level mapping of 
UniProt, Pfam, and other databases with PDB entries.

- [Pfam](@ref Module-Pfam): This module uses the previous modules to work with Pfam MSAs. 
It also has useful parameter optimization functions to be used with Pfam alignments.

- [Utils](@ref API-Utils): MIToS has also a Utils module with common utils functions and 
types used in different modules of this package.

## Citation  

If you use MIToS, please cite:

Diego J. Zea, Diego Anfossi, Morten Nielsen, Cristina Marino-Buslje; MIToS.jl: mutual information tools for protein sequence analysis in the Julia language, Bioinformatics, Volume 33, Issue 4, 15 February 2017, Pages 564–565, [https://doi.org/10.1093/bioinformatics/btw646](https://doi.org/10.1093/bioinformatics/btw646)  

## Older MIToS versions

You can change the MIToS version of the documentation at the bottom left of this site—the 
older version available is MIToS 2.0. If you are using MIToS v1 in a version of Julia 
pre-1.0, please read [this older documentation](https://diegozea.github.io/mitosghpage-legacy/) instead.  


![Leloir Institute Foundation](./assets/FIL_logo.png)  

Structural Bioinformatics Unit, Leloir Institute Foundation.
Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
