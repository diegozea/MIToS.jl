
# Scripts

MIToS implements several useful scripts to **command line execution (without requiring Julia coding)**. All this scripts are located in the `scripts` folder of the MIToS directory. You can copy them to your working directory, used the path to their folder or put them in the path (look into the **Installation** section of this manual).

## Contents

- [Corrected Mutual Information](#Corrected-Mutual-Information)
    - [Buslje09.jl](#Buslje09.jl)
    - [BLMI.jl](#BLMI.jl)
- [Protein Structure](#Protein-Structure) 
    - [DownloadPDB.jl](#DownloadPDB.jl)
    - [Distances.jl](#Distances.jl)
- [MSA](#MSA)
    - [MSADescription.jl](#MSADescription.jl)
    - [PercentIdentity.jl](#PercentIdentity.jl)
    - [AlignedColumns.jl](#AlignedColumns.jl)
    - [SplitStockholm.jl](#SplitStockholm.jl)

## Corrected Mutual Information

<a href="#"><i class="fa fa-arrow-up"></i></a>

### Buslje09.jl


```julia
;Buslje09.jl -h
```

    usage: Buslje09.jl [-l] [-o OUTPUT] [-p PARALLEL] [-f FORMAT]
                       [-L LAMBDA] [-c] [-i THRESHOLD] [-g MAXGAP] [-a]
                       [-s SAMPLES] [-G] [-F] [--version] [-h] [FILE]
    
    This takes a MSA file as input. Calculates and saves on *.busjle09.csv
    a Z score and a corrected MI/MIp as described on:
    Buslje, C. M., Santos, J., Delfino, J. M., & Nielsen, M. (2009).
    Correction for phylogeny, small number of observations and data
    redundancy improves the identification of coevolving amino acid pairs
    using mutual information. Bioinformatics, 25(9), 1125-1131.
    
    positional arguments:
      FILE                  File name. If it is not used, the script reads
                            from STDIN.
    
    optional arguments:
      -l, --list            The input is a list of file names. If -p is
                            used, files will be processed in parallel.
      -o, --output OUTPUT   Name of the output file. Output will be gzip
                            if the extension is ".gz". If it starts with a
                            dot, the name is used as a suffix or extension
                            of the input filename. If it ends with a dot,
                            is used as a prefix. If the output name starts
                            and ends with dots, it's used as an interfix
                            before the extension.If a single file is used
                            and there is not a file name (STDIN), the
                            output will be print into STDOUT, unless a
                            output filename is used. You can use "STDOUT"
                            to force print into STDOUT. STDOUT can not be
                            use with --list. (default: ".busjle09.csv")
      -p, --parallel PARALLEL
                            Number of worker processes. (type: Int64,
                            default: 1)
      -f, --format FORMAT   Format of the MSA: Stockholm, Raw or FASTA
                            (default: "Stockholm")
      -L, --lambda LAMBDA   Low count value (type: Float64, default: 0.05)
      -c, --clustering      Sequence clustering (Hobohm I)
      -i, --threshold THRESHOLD
                            Percent identity threshold for sequence
                            clustering (Hobohm I) (type: Float64, default:
                            62.0)
      -g, --maxgap MAXGAP   Maximum fraction of gaps in positions included
                            in calculation (type: Float64, default: 0.5)
      -a, --apc             Use APC correction (MIp)
      -s, --samples SAMPLES
                            Number of samples for Z-score (type: Int64,
                            default: 100)
      -G, --usegap          Use gaps on statistics
      -F, --fixedgaps       Fix gaps positions for the random samples
      --version             show version information and exit
      -h, --help            show this help message and exit
    
    
    MIToS 1.1.0+
    
    Bioinformatics Unit
    Leloir Institute Foundation
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    


<a href="#"><i class="fa fa-arrow-up"></i></a>

### BLMI.jl


```julia
;BLMI.jl -h
```

    usage: BLMI.jl [-l] [-o OUTPUT] [-p PARALLEL] [-f FORMAT] [-b BETA]
                   [-i THRESHOLD] [-g MAXGAP] [-a] [-s SAMPLES] [-F]
                   [--version] [-h] [FILE]
    
    This takes a MSA file as input. Calculates and saves on *.BLMI.csv a Z
    score and a corrected MI/MIp. The script uses BLOSUM62 based pseudo
    frequencies and sequences clustering (Hobohm I).
    
    positional arguments:
      FILE                  File name. If it is not used, the script reads
                            from STDIN.
    
    optional arguments:
      -l, --list            The input is a list of file names. If -p is
                            used, files will be processed in parallel.
      -o, --output OUTPUT   Name of the output file. Output will be gzip
                            if the extension is ".gz". If it starts with a
                            dot, the name is used as a suffix or extension
                            of the input filename. If it ends with a dot,
                            is used as a prefix. If the output name starts
                            and ends with dots, it's used as an interfix
                            before the extension.If a single file is used
                            and there is not a file name (STDIN), the
                            output will be print into STDOUT, unless a
                            output filename is used. You can use "STDOUT"
                            to force print into STDOUT. STDOUT can not be
                            use with --list. (default: ".BLMI.csv")
      -p, --parallel PARALLEL
                            Number of worker processes. (type: Int64,
                            default: 1)
      -f, --format FORMAT   Format of the MSA: Stockholm, Raw or FASTA
                            (default: "Stockholm")
      -b, --beta BETA       β for BLOSUM62 pseudo frequencies (type:
                            Float64, default: 8.512)
      -i, --threshold THRESHOLD
                            Percent identity threshold for sequence
                            clustering (Hobohm I) (type: Float64, default:
                            62.0)
      -g, --maxgap MAXGAP   Maximum fraction of gaps in positions included
                            in calculation (type: Float64, default: 0.5)
      -a, --apc             Use APC correction (MIp)
      -s, --samples SAMPLES
                            Number of samples for Z-score (type: Int64,
                            default: 50)
      -F, --fixedgaps       Fix gaps positions for the random samples
      --version             show version information and exit
      -h, --help            show this help message and exit
    
    
    MIToS 1.1.0+
    
    Bioinformatics Unit
    Leloir Institute Foundation
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    


<a href="#"><i class="fa fa-arrow-up"></i></a>

## Protein Structure

<a href="#"><i class="fa fa-arrow-up"></i></a>

### DownloadPDB.jl


```julia
;DownloadPDB.jl -h
```

    usage: DownloadPDB.jl [-c CODE] [-l LIST] [-f FORMAT] [--version] [-h]
    
    Download gzipped files from PDB.
    
    optional arguments:
      -c, --code CODE      PDB code
      -l, --list LIST      File with a list of PDB codes (one per line)
      -f, --format FORMAT  Format. It should be PDBFile (pdb) or PDBML
                           (xml) (default: "PDBML")
      --version            show version information and exit
      -h, --help           show this help message and exit
    
    
    MIToS 1.1.0+
    
    Bioinformatics Unit
    Leloir Institute Foundation
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    


<a href="#"><i class="fa fa-arrow-up"></i></a>

### Distances.jl


```julia
;Distances.jl -h
```

    usage: Distances.jl [-l] [-o OUTPUT] [-p PARALLEL] [-d DISTANCE]
                        [-f FORMAT] [-m MODEL] [-c CHAIN] [-g GROUP] [-i]
                        [--version] [-h] [FILE]
    
    Calculates residues distance and writes them into a *.distances.csv.gz
    gzipped file.
    
    positional arguments:
      FILE                  File name. If it is not used, the script reads
                            from STDIN.
    
    optional arguments:
      -l, --list            The input is a list of file names. If -p is
                            used, files will be processed in parallel.
      -o, --output OUTPUT   Name of the output file. Output will be gzip
                            if the extension is ".gz". If it starts with a
                            dot, the name is used as a suffix or extension
                            of the input filename. If it ends with a dot,
                            is used as a prefix. If the output name starts
                            and ends with dots, it's used as an interfix
                            before the extension.If a single file is used
                            and there is not a file name (STDIN), the
                            output will be print into STDOUT, unless a
                            output filename is used. You can use "STDOUT"
                            to force print into STDOUT. STDOUT can not be
                            use with --list. (default:
                            ".distances.csv.gz")
      -p, --parallel PARALLEL
                            Number of worker processes. (type: Int64,
                            default: 1)
      -d, --distance DISTANCE
                            The distance to be calculated, options: All,
                            Heavy, CA, CB (default: "All")
      -f, --format FORMAT   Format of the PDB file: It should be PDBFile
                            or PDBML (default: "PDBFile")
      -m, --model MODEL     The model to be used, * for all (default: "1")
      -c, --chain CHAIN     The chain to be used, * for all (default: "*")
      -g, --group GROUP     Group of atoms to be used, should be ATOM,
                            HETATM or * for all (default: "*")
      -i, --inter           Calculate inter chain distances
      --version             show version information and exit
      -h, --help            show this help message and exit
    
    
    MIToS 1.1.0+
    
    Bioinformatics Unit
    Leloir Institute Foundation
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    


<a href="#"><i class="fa fa-arrow-up"></i></a>

## MSA

<a href="#"><i class="fa fa-arrow-up"></i></a>

### MSADescription.jl


```julia
;MSADescription.jl -h
```

    usage: MSADescription.jl [-l] [-o OUTPUT] [-p PARALLEL] [-f FORMAT]
                            [-e] [--version] [-h] [FILE]
    
    Creates an *.description.csv from a Stockholm file with: the number of
    columns, sequences, clusters after Hobohm clustering at 62% identity
    and mean percent identity. Also the mean, standard deviation and
    quantiles of: sequence coverage of the MSA, gap percentage.
    
    positional arguments:
      FILE                  File name. If it is not used, the script reads
                            from STDIN.
    
    optional arguments:
      -l, --list            The input is a list of file names. If -p is
                            used, files will be processed in parallel.
      -o, --output OUTPUT   Name of the output file. Output will be gzip
                            if the extension is ".gz". If it starts with a
                            dot, the name is used as a suffix or extension
                            of the input filename. If it ends with a dot,
                            is used as a prefix. If the output name starts
                            and ends with dots, it's used as an interfix
                            before the extension.If a single file is used
                            and there is not a file name (STDIN), the
                            output will be print into STDOUT, unless a
                            output filename is used. You can use "STDOUT"
                            to force print into STDOUT. STDOUT can not be
                            use with --list. (default: ".description.csv")
      -p, --parallel PARALLEL
                            Number of worker processes. (type: Int64,
                            default: 1)
      -f, --format FORMAT   Format of the MSA: Stockholm, Raw or FASTA
                            (default: "Stockholm")
      -e, --exact           If it's true, the mean percent identity is
                            exact (using all the pairwise comparisons).
      --version             show version information and exit
      -h, --help            show this help message and exit
    
    
    MIToS 1.1.0+
    
    Bioinformatics Unit
    Leloir Institute Foundation
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    


<a href="#"><i class="fa fa-arrow-up"></i></a>

### PercentIdentity.jl


```julia
;PercentIdentity.jl -h
```

    usage: PercentIdentity.jl [-l] [-o OUTPUT] [-p PARALLEL] [-f FORMAT]
                            [-s] [--version] [-h] [FILE]
    
    Calculates the percentage identity between all the sequences of an MSA
    and creates an *.pidstats.csv file with: The number of columns and
    sequences. The mean, standard deviation, median, minimum and maximum
    values and first and third quantiles of the percentage identity. It
    could also create and pidlist.csv file with the percentage identity
    for each pairwise comparison.
    
    positional arguments:
      FILE                  File name. If it is not used, the script reads
                            from STDIN.
    
    optional arguments:
      -l, --list            The input is a list of file names. If -p is
                            used, files will be processed in parallel.
      -o, --output OUTPUT   Name of the output file. Output will be gzip
                            if the extension is ".gz". If it starts with a
                            dot, the name is used as a suffix or extension
                            of the input filename. If it ends with a dot,
                            is used as a prefix. If the output name starts
                            and ends with dots, it's used as an interfix
                            before the extension.If a single file is used
                            and there is not a file name (STDIN), the
                            output will be print into STDOUT, unless a
                            output filename is used. You can use "STDOUT"
                            to force print into STDOUT. STDOUT can not be
                            use with --list. (default: ".pidstats.csv")
      -p, --parallel PARALLEL
                            Number of worker processes. (type: Int64,
                            default: 1)
      -f, --format FORMAT   Format of the MSA: Stockholm, Raw or FASTA
                            (default: "Stockholm")
      -s, --savelist        Create and pidlist.csv file with the
                            percentage identity for each pairwise
                            comparison.
      --version             show version information and exit
      -h, --help            show this help message and exit
    
    
    MIToS 1.1.0+
    
    Bioinformatics Unit
    Leloir Institute Foundation
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    


<a href="#"><i class="fa fa-arrow-up"></i></a>

### AlignedColumns.jl


```julia
;AlignedColumns.jl -h
```

    usage: AlignedColumns.jl [-l] [-o OUTPUT] [-p PARALLEL] [--version]
                            [-h] [FILE]
    
    Creates a file in Stockholm format with the aligned columns from a
    Pfam Stockholm file. Insertions are deleted, as they are unaligned in
    a proﬁle HMM. The output file *.aligned.* contains UniProt residue
    numbers and original column numbers in its annotations.
    
    positional arguments:
      FILE                  File name. If it is not used, the script reads
                            from STDIN.
    
    optional arguments:
      -l, --list            The input is a list of file names. If -p is
                            used, files will be processed in parallel.
      -o, --output OUTPUT   Name of the output file. Output will be gzip
                            if the extension is ".gz". If it starts with a
                            dot, the name is used as a suffix or extension
                            of the input filename. If it ends with a dot,
                            is used as a prefix. If the output name starts
                            and ends with dots, it's used as an interfix
                            before the extension.If a single file is used
                            and there is not a file name (STDIN), the
                            output will be print into STDOUT, unless a
                            output filename is used. You can use "STDOUT"
                            to force print into STDOUT. STDOUT can not be
                            use with --list. (default: ".aligned.")
      -p, --parallel PARALLEL
                            Number of worker processes. (type: Int64,
                            default: 1)
      --version             show version information and exit
      -h, --help            show this help message and exit
    
    
    MIToS 1.1.0+
    
    Bioinformatics Unit
    Leloir Institute Foundation
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    


<a href="#"><i class="fa fa-arrow-up"></i></a>

### SplitStockholm.jl


```julia
;SplitStockholm.jl -h
```

    usage: SplitStockholm.jl [-p PATH] [--version] [-h] file
    
    Splits a file with multiple sequence alignments in Stockholm format,
    creating one compressed file per MSA in Stockholm format:
    accessionumber.gz
    
    positional arguments:
      file             Input file
    
    optional arguments:
      -p, --path PATH  Path for the output files [default: execution
                       directory] (default: "")
      --version        show version information and exit
      -h, --help       show this help message and exit
    
    
    MIToS 1.1.0+
    
    Bioinformatics Unit
    Leloir Institute Foundation
    Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina
    


<a href="#"><i class="fa fa-arrow-up"></i></a>
