var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "(Image: MIToS.jl)  MIToS is an environment for Mutual Information (MI) analysis and implements several useful tools for Multiple Sequence Alignments (MSAs) and PDB structures management in the Julia language.This is the documentation for MIToS 2.0 in Julia 0.5. If you are using MIToS 1.0 in Julia 0.4, please read this documentation  instead.  "
},

{
    "location": "index.html#Modules-1",
    "page": "Home",
    "title": "Modules",
    "category": "section",
    "text": "MIToS tools are separated on different modules, related to different tasks.MSA : This module defines multiple functions and types for dealing with MSAs and their annotations. It also includes facilities for sequence clustering.PDB : This module defines types and methods to work with protein structures from PDB.SIFTS : This module allows access to SIFTS residue-level mapping of UniProt, Pfam and other databases with PDB entries.Information : This module defines residue contingency tables and methods on them to estimate information measure from MSAs. It includes functions to estimate corrected mutual information (ZMIp, ZBLMIp) between MSA columns.  Pfam : This module use the previous modules to work with Pfam MSAs. It also has useful parameter optimization functions to be used with Pfam alignments.  Utils : MIToS has also an Utils module with common utils functions and types used in this package.  "
},

{
    "location": "index.html#Citation-1",
    "page": "Home",
    "title": "Citation",
    "category": "section",
    "text": "If you use MIToS, please cite:Diego J. Zea, Diego Anfossi, Morten Nielsen, Cristina Marino-Buslje; MIToS.jl: mutual information tools for protein sequence analysis in the Julia language, Bioinformatics, Volume 33, Issue 4, 15 February 2017, Pages 564–565, https://doi.org/10.1093/bioinformatics/btw646  (Image: Leloir Institute Foundation)  Structural Bioinformatics Unit, Leloir Institute Foundation. Av. Patricias Argentinas 435, CP C1405BWE, Buenos Aires, Argentina"
},

{
    "location": "Installation.html#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "Installation.html#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "First you need Julia 0.5.2(Image: ) installed. A new version of MIToS with Julia 0.6 support is under development, we expect to release it soon. The MIToS' stable version can be installed by typing on the Julia REPL:  Pkg.add(\"MIToS\")If everything goes well with the installation, MIToS will be loaded without errors by typing:  using MIToSYou can optionally do an exhaustive test of your installed version of MIToS with Pkg.test (it takes few seconds):  Pkg.test(\"MIToS\")note: Note\nWays to run Julia  Option Description\nJulia REPL(Image: ) Built-in Julia command line. Start an Julia interactive session (REPL) by double-clicking the Julia executable or running julia from the system command line.\nJuliaBox(Image: ) You can try Julia from your web browser. No installation is required.\nIJulia(Image: ) Jupyter/IPython notebook for Julia. It was used for generating the this documentation.\nJuno(Image: ) Integrated Development Environment (IDE)."
},

{
    "location": "Installation.html#Plots-installation-1",
    "page": "Installation",
    "title": "Plots installation",
    "category": "section",
    "text": "Julia plotting capabilities are available through external packages. MIToS make use of  RecipesBase to define plot recipes, which can be plotted using  Plots(Image: ) and different  backends. You need to install Plots(Image: )  to plot MIToS objects:  Pkg.add(\"Plots\")And you also need to install at least one of the following backends:  Pkg.add(\"PyPlot\")\nPkg.add(\"GR\") # Fast\nPkg.add(\"PlotlyJS\") # InteractiveYou need to load Plots in order to use the plot function. There is more information about it in the Plots documentation(Image: ).  using PlotsTo generate graph (network), arc and chord (circo) plots, you also need to install and load PlotRecipes(Image: ).  Pkg.add(\"PlotRecipes\")\n\nusing PlotRecipes"
},

{
    "location": "Installation.html#Scripts-location-1",
    "page": "Installation",
    "title": "Scripts location",
    "category": "section",
    "text": "The MIToS’ scripts are located in the MIToS/scripts folder and can be runned from your system command line. It’s possible to ask Julia for the location of the installed package using Pkg.dir  joinpath(Pkg.dir(\"MIToS\"), \"scripts\")You might want to add this folder into your PATH to easily access MIToS’ scripts.  "
},

{
    "location": "Installation.html#How-to-add-the-script-folder-to-PATH-in-Bash?-1",
    "page": "Installation",
    "title": "How to add the script folder to PATH in Bash?",
    "category": "section",
    "text": "You can do it by adding the path of the MIToS script folder into the ~/.bashrc file:open(joinpath(homedir(), \".bashrc\"), \"r+\") do fh\n    path_to_scripts = joinpath(Pkg.dir(\"MIToS\"), \"scripts\")\n    if all(line -> !contains(line, path_to_scripts), eachline(fh))\n        println(fh, \"export PATH=\\\"\\$PATH:\", path_to_scripts, \"\\\"\")\n    end\nend"
},

{
    "location": "Example.html#",
    "page": "Example",
    "title": "Example",
    "category": "page",
    "text": ""
},

{
    "location": "Example.html#Example-1",
    "page": "Example",
    "title": "Example",
    "category": "section",
    "text": "In this simple demonstration, you will see how to calculate ZBLMIp (Z score of the corrected MIp using BLOSUM62 pseudo frequencies) for a Pfam(Image: ) MSA from the Julia REPL or using a MIToS script in the system command line.  "
},

{
    "location": "Example.html#juliarepl-1",
    "page": "Example",
    "title": "MIToS in the Julia REPL",
    "category": "section",
    "text": "If you load the Pfam module from MIToS, you will get access to a set of functions that work with Pfam MSAs. In this case, we are going to use it for download a Stockholm(Image: ) MSA from the Pfam website and read it into Julia.  using Plots\npyplot() # Just to avoid warnings in the outputusing MIToS.Pfam\npfam_file = downloadpfam(\"PF10660\")\nmsa = read(pfam_file, Stockholm, generatemapping=true, useidcoordinates=true)note: Note\nGeneration of sequence and column mappings   The keyword argument generatemapping of read allows to generate sequence and column mappings for the MSA. Column mapping is the map between of each column on the MSA object and the column number in the file. Sequence mappings will use the start and end coordinates in the sequence ids for enumerate each residue in the sequence if useidcoordinates is true.  You can plot this MSA and other MIToS’ objects using the Plots(Image: ) package. The installation of Plots is described in the Installation section of this site:using Plots\npyplot() # Best backend for MSA plots\nplot(msa)\npng(\"msa.png\") # hide\nnothing # hide(Image: )  The Information module of MIToS has functions to calculate measures from the Information Theory(Image: ), such as Entropy and Mutual Information (MI), on a MSA. In this example, we will estimate covariation between columns of the MSA with a corrected MI that use the BLOSUM62 matrix for calculate pseudo frequencies (BLMI).  using MIToS.Information\nZBLMIp, BLMIp = BLMI(msa)\nZBLMIp # shows ZBLMIp scoresOnce the Plots package is installed and loaded, you can use its capabilities to visualize this results:heatmap(ZBLMIp, yflip=true, c=:grays)\npng(\"blmi.png\") # hide\nnothing # hide(Image: )  rm(pfam_file) # clean up"
},

{
    "location": "Example.html#commandline-1",
    "page": "Example",
    "title": "MIToS in system command line",
    "category": "section",
    "text": "Calculate ZBLMIp on the system shell is easy using the MIToS script called BLMI.jl. This script reads a MSA file, and writes a file with the same base name of the input but with the .BLMI.csv extension.  BLMI.jl PF14972.stockholm.gz"
},

{
    "location": "MSA.html#",
    "page": "MSA",
    "title": "MSA",
    "category": "page",
    "text": ""
},

{
    "location": "MSA.html#Module-MSA-1",
    "page": "MSA",
    "title": "MSA",
    "category": "section",
    "text": "The MSA module of MIToS has utilities for working with Multiple Sequence Alignments of protein Sequences (MSA).  using MIToS.MSA # to load the MSA module"
},

{
    "location": "MSA.html#Features-1",
    "page": "MSA",
    "title": "Features",
    "category": "section",
    "text": "Read and write MSAs in Stockholm, FASTA or Raw format.\nHandle MSA annotations.\nEdit the MSA, e.g. delete columns or sequences, change sequence order, shuffling...\nKeep track of positions and annotations after modifications on the MSA.\nDescribe an MSA, e.g. mean percent identity, sequence coverage, gap percentage...\nSequence clustering with Hobohm I."
},

{
    "location": "MSA.html#Contents-1",
    "page": "MSA",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"MSA.md\"]\nDepth = 4"
},

{
    "location": "MSA.html#MSA-IO-1",
    "page": "MSA",
    "title": "MSA IO",
    "category": "section",
    "text": ""
},

{
    "location": "MSA.html#Reading-MSA-files-1",
    "page": "MSA",
    "title": "Reading MSA files",
    "category": "section",
    "text": "The main function for reading MSA files in MIToS is read and it is defined in the Utils module. This function takes a filename/path as a first argument followed by other arguments. It opens the file and uses the arguments to call the parse function. read decides how to open the file, using the prefixes and suffixes of the file name, while parse does the actual parsing of the file. You can read gzipped files if they have the .gz extension and also urls pointing to a web file.   The second argument of read and parse is the file Format. The supported MSA formats at the moment are Stockholm, FASTA and Raw.   For example, reading with MIToS the full Stockholm MSA of the family PF07388 using the Pfam RESTful interface will be:  using MIToS.MSA\n\nread(\"http://pfam.xfam.org/family/PF07388/alignment/full\", Stockholm)The third (and optional) argument of read and parse is the output MSA type:  Matrix{Residue} : It contains the aligned sequences.  \nMultipleSequenceAlignment : It contains the aligned sequences and theirnames/identifiers.  AnnotatedMultipleSequenceAlignment : Is the richest MSA format of MIToS. It's thedefault. It includes the aligned sequences, their names and the MSA annotations.  Example of Matrix{Residue} output using a Stockholm file as input:read(\"http://pfam.xfam.org/family/PF07388/alignment/full\", Stockholm, Matrix{Residue})Given that read calls parse, you should look into the documentation of the last one to know the available keyword arguments. The optional keyword arguments of those functions are:  generatemapping : If generatemapping is true (default: false), sequences andcolumns mappings are generated and saved in the MSA annotations. The default is false to not overwrite mappings by mistake when you read an annotated MSA file saved with MIToS.  useidcoordinates : If useidcoordinates is true (default: false) and the nameshave the form seqname/start-end, MIToS uses this coordinates to generate sequence mappings. This is safe and useful with unmodified Pfam MSAs. Do not use it when reading an MSA saved with MIToS. MIToS deletes unaligned insert columns, therefore disrupts sequences that have them.  deletefullgaps : Given that lowercase characters and dots are converted to gaps,unaligned insert columns in the MSA (derived from a HMM profile) are converted into full gap columns. deletefullgaps is true by default, deleting full gaps columns and therefore insert columns.  note: Note\nIf you want to keep the insert columns...  Use the keyword argument keepinserts to true in read/parse. This only works with an AnnotatedMultipleSequenceAlignment output. A column annotation (\"Aligned\") is stored in the annotations, where insert columns are marked with 0 and aligned columns with 1.  When read returns an AnnotatedMultipleSequenceAlignment, it uses the MSA Annotations to keep track of performed modifications. To access this notes, use printmodifications:  msa = read(\"http://pfam.xfam.org/family/PF01565/alignment/full\", Stockholm)\n\nprintmodifications(msa)"
},

{
    "location": "MSA.html#Writing-MSA-files-1",
    "page": "MSA",
    "title": "Writing MSA files",
    "category": "section",
    "text": "Julia REPL shows MSAs as Matrices. If you want to print them in another format, you should use the print function with an MSA object as first argument and the Format FASTA, Stockholm or Raw as second argument.  using MIToS.MSA\n\nmsa = read(\"http://pfam.xfam.org/family/PF16996/alignment/full\", Stockholm) # reads a Stockholm MSA file\n\nprint(msa, FASTA) # prints msa in FASTA formatTo save an MSA object to a file, use the write function. This function takes a filename as a first argument. If the filename ends with .gz, the output will be a compressed (gzipped) file. The next two arguments of write are passed to print, so write behaves as print.  write(\"msa.gz\", msa, FASTA) # writes msa in FASTA format in a gzipped file"
},

{
    "location": "MSA.html#MSA-Annotations-1",
    "page": "MSA",
    "title": "MSA Annotations",
    "category": "section",
    "text": "MSA annotations are based on the Stockholm format mark-ups. There are four types of annotations stored as dictionaries. All the annotations have a feature name as part of the key, which should be a single \"word\" (without spaces) and less than 50 characters long.  File annotations : The annotations can contain either file or MSA information. Theyhave feature names as keys and the values are strings (free text). Lines starting with #=GF in Stockholm format.  Column annotations : They have feature names as keys and strings with exactly 1 charper column as values. Lines starting with #=GC in Stockholm format.  Sequence annotations : The keys are tuples with the sequence name and the featurename. The values are free text (strings). Lines starting with #=GS in Stockholm format.  Residue annotations : The keys are tuples with the sequence name and the featurename. The values are strings with exactly 1 char per column/residues. #=GR lines in Stockholm format.  Julia REPL shows the Annotations type as they are represented in the Stockholm format(Image: ). You can get the Annotations inside an annotated MSA or sequence using the annotations function.  annotations(msa)Particular annotations can be accessed using the functions getannot.... This functions take the MSA/sequence as first argument and the feature name of the desired annotation as the last. In the case of getannotsequence and getannotresidue, the second argument should be the sequence name.  getannotsequence(msa, \"J0UVX5_STREE/3-39\", \"AC\") # (\"J0UVX5_STREE/3-39\", \"AC\") is the key in the dictionaryIf you want to add new annotations, you should use the setannot…! functions. This functions have the same arguments that getannot... functions but with an extra argument to indicate the new annotation.  setannotsequence!(msa, \"J0UVX5_STREE/3-39\", \"New_Feature_Name\", \"New_Annotation\")A getannot... function without the key (last arguments), returns the particular annotation dictionary. As you can see, the new sequence annotation is now part of our MSA annotations.  getannotsequence(msa)"
},

{
    "location": "MSA.html#Editing-your-MSA-1",
    "page": "MSA",
    "title": "Editing your MSA",
    "category": "section",
    "text": "MIToS offers functions to edit your MSA. Given that this functions modify the msa, their names end with a bang !, following the Julia convention. Some of these functions have an annotate keyword argument (in general it is true by default) to indicate if the modification should be recorded in the MSA/sequence annotations.  One common task is to delete sequences or columns of the MSA. This could be done using the functions filtersequences! and filtercolumns!. This functions take the MSA or sequence (if it's possible) as first argument and a BitVector or Vector{Bool} mask as second argument. It deletes all the sequences or columns where the mask is false. This functions are also defined for Annotations, this allows to automatically update (modify) the annotations (and therefore, sequence and column mappings) in the MSA.  This two deleting operations are used in the second and third mutating functions of the following list:  setreference! : Sets one of the sequences as the first sequence of the MSA (query orreference sequence).  adjustreference! : Deletes columns with gaps in the first sequence of the MSA(reference).  gapstrip! : This function first calls adjustreference!, then deletes sequences withlow (user defined) MSA coverage and finally, columns with user defined % of gaps.  Also, there are several available funtions shuffle_…!. These functions are useful to generate random alignments. The Information module of MIToS uses them to calculate the Z scores of MI values.  "
},

{
    "location": "MSA.html#Example:-Deleting-sequences-1",
    "page": "MSA",
    "title": "Example: Deleting sequences",
    "category": "section",
    "text": "For example, if you want to keep only the proteins from Actinobacteria you can delete all the sequences that don't have _9ACTN in their UniProt entry names:  using MIToS.MSA\n\nmsa = read(\"http://pfam.xfam.org/family/PF07388/alignment/full\", Stockholm)\n\nsequencenames(msa) # the function sequencenames returns the sequence names in the MSAmask = map(x -> ismatch(r\"_9ACTN\", x), sequencenames(msa)) # an element of mask is true if \"_9ACTN\" is in the namefiltersequences!(msa, mask) # deletes all the sequences where mask is false\n\nsequencenames(msa)"
},

{
    "location": "MSA.html#Example:-Exporting-a-MSA-for-freecontact-(part-I)-1",
    "page": "MSA",
    "title": "Example: Exporting a MSA for freecontact (part I)",
    "category": "section",
    "text": "The most simple input for the command line tool freecontact(Image: ) (if you don't want to set --mincontsep) is a Raw MSA file with a reference sequence without insertions or gaps. This is easy to get with MIToS using read (deletes the insert columns), setreference! (to choose a reference), adjustreference! (to delete columns with gaps in the reference) and write (to save it in Raw format) functions.  using MIToS.MSA\nmsa = read(\"http://pfam.xfam.org/family/PF02476/alignment/full\", Stockholm)\nmaxcoverage, indice = findmax(coverage(msa)) # chooses the sequence with more coverage of the MSA\nsetreference!(msa, indice)\nadjustreference!(msa)\nwrite(\"tofreecontact.msa\", msa, Raw)\nprint(readstring(\"tofreecontact.msa\")) # It displays the contents of the output file"
},

{
    "location": "MSA.html#Column-and-sequence-mappings-1",
    "page": "MSA",
    "title": "Column and sequence mappings",
    "category": "section",
    "text": "Inserts in a Stockholm MSA allow to access the full fragment of the aligned sequences. Using this, combined with the sequence names with coordinates used in Pfam, you can know what is the UniProt residue number of each residue in the MSA.   \"PROT_SPECI/3-15 .....insertALIGNED\"\n#                     3456789111111\n#                            012345MIToS read and parse functions deletes the insert columns, but they do the mapping of each residue to its residue number before deleting insert columns if generatemapping is true. If you don't set useidcoordinates to true, the residue first i residue will be 1 instead of 3 in the previous example.  using MIToS.MSA\n\nmsa = parse(\"PROT_SPECI/3-15 .....insertALIGNED\", Stockholm, generatemapping=true, useidcoordinates=true)MIToS also keeps the column number of the input MSA and its total number of columns. All this data is stored in the MSA annotations using the SeqMap, ColMap and NCol feature names.  annotations(msa)To have an easy access to mapping data, MIToS provides the getsequencemapping and getcolumnmapping functions.  getsequencemapping(msa, \"PROT_SPECI/3-15\")getcolumnmapping(msa)"
},

{
    "location": "MSA.html#Example:-Exporting-a-MSA-for-freecontact-(part-II)-1",
    "page": "MSA",
    "title": "Example: Exporting a MSA for freecontact (part II)",
    "category": "section",
    "text": "If we want to use the --mincontsep argument of freecontact to calculate scores between distant residues, we will need to add a header to the MSA. This header should contains the residue number of the first residue of the sequence and the full fragment of that sequence (with the inserts). This data is used by FreeContact to calculate the residue number of each residue in the reference sequence.   We are going to use MIToS mapping data to create this header, so we read the MSA with generatemapping and useidcoordinates setted to true.  using MIToS.MSA\n\nmsa = read( \"http://pfam.xfam.org/family/PF02476/alignment/full\", Stockholm,\n            generatemapping=true, useidcoordinates=true)Here, we are going to choose the sequence with more coverage of the MSA as our reference sequence.  maxcoverage, indice = findmax(coverage(msa))\nsetreference!(msa, indice)\nadjustreference!(msa)MIToS deletes the residues in insert columns, so we are going to use the sequence mapping to generate the whole fragment of the reference sequence (filling the missing regions with 'x').  seqmap = getsequencemapping(msa, 1) # seqmap will be a vector with the residue numbers of the first sequence (reference)\n\nseq = collect( stringsequence(msa, 1) ) # seq will be a Vector of Chars with the reference sequence\n\nsequence = map(seqmap[1]:seqmap[end]) do seqpos # for each position in the whole fragment\n    if seqpos in seqmap                         # if that position is in the MSA\n        shift!(seq)                             # the residue is taken from seq\n    else                                        # otherwise\n        'x'                                     # 'x' is included\n    end\nend\n\nsequence = join(sequence) # join the Chars on the Vector to create a stringOnce we have the whole fragment of the sequence, we create the file and write the header in the required format (as in the man page of freecontact).  open(\"tofreecontact.msa\", \"w\") do fh\n\n    println(fh, \"# querystart=\", seqmap[1])\n\n    println(fh, \"# query=\", sequence )\n\nendAs last (optional) argument, write takes the mode in which is opened the file. We use \"a\" here to append the MSA to the header.  write(\"tofreecontact.msa\", msa, Raw, \"a\")print(readstring(\"tofreecontact.msa\")) # It displays the contents of the output file"
},

{
    "location": "MSA.html#Get-sequences-from-a-MSA-1",
    "page": "MSA",
    "title": "Get sequences from a MSA",
    "category": "section",
    "text": "It's possible to index the MSA as any other matrix to get an aligned sequence. This will be return a Array of Residues without annotations but keeping names/identifiers.  using MIToS.MSA\n\nmsa = read( \"http://pfam.xfam.org/family/PF16996/alignment/full\", Stockholm,\n            generatemapping=true, useidcoordinates=true)msa[2,:] # second sequence of the MSA, it keeps column namesmsa[2:2,:] # Using the range 2:2 to select the second sequence, keeping the sequence nameIf you want to obtain the aligned sequence with its name and annotations (and therefore sequence and column mappings), you should use the function getsequence. This function returns an AlignedSequence with the sequence name from a MultipleSequenceAlignment or an AnnotatedAlignedSequence, that also contains annotations, from an AnnotatedMultipleSequenceAlignment.  secondsequence = getsequence(msa, 2)annotations(secondsequence)Use stringsequence if you want to get the sequence as a string.  stringsequence(msa, 2)Given that matrices are stored columnwise in Julia, you will find useful the getresiduesequences function when you need to heavily operate over sequences.  getresiduesequences(msa)"
},

{
    "location": "MSA.html#Describing-your-MSA-1",
    "page": "MSA",
    "title": "Describing your MSA",
    "category": "section",
    "text": "The MSA module has a number of functions to gain insight about your MSA. Using MIToS.MSA, one can easily ask for...  The number of columns and sequences with the ncolumns and nsequences functions.  \nThe fraction of columns with residues (coverage) for each sequence making use of thecoverage method.  The fraction or percentage of gaps/residues using with the functions gapfraction,residuefraction and columngapfraction.  The percentage of identity (PID) between each sequence of the MSA or its mean valuewith percentidentity and meanpercentidentity.  The percentage identity between two aligned sequences it's a common measure of sequence similarity and it's used by the hobohmI method to estimate and reduce MSA redundancy. MIToS functions to calculate percent identity don't align the sequences, they need sequences already aligned. Full gaps columns don't count to the align length.  using MIToS.MSA\n\nmsa = permutedims(\n        hcat(   res\"--GGG-\",      # res\"...\" uses the @res_str macro to create a (column) Vector{Residue}\n                res\"---GGG\" ), (2,1))\n#        identities 000110 sum 2\n#  aligned residues 001111 sum 4percentidentity(msa[1,:], msa[2,:]) # 2 / 4To quickly calculate if the percentage of identity is greater than a determined value, use that threshold as third argument. percentidentity(seqa, seqb, pid) is a lot more faster than percentidentity(seqa, seqb) >= pid.  percentidentity(msa[1,:], msa[2,:], 62) # 50% >= 62%"
},

{
    "location": "MSA.html#Example:-Plotting-gap-percentage-per-column-and-coverage-per-sequence-1",
    "page": "MSA",
    "title": "Example: Plotting gap percentage per column and coverage per sequence",
    "category": "section",
    "text": "The gapfraction and coverage functions return a vector of number between 0.0 and 1.0 (fraction of...). Sometime it's useful to plot this data to quickly understand the MSA structure. In this example, we are going to use the Plots(Image: ) package for plotting, with a PyPlot(Image: ) backend, but you are free to use any of the Julia plotting libraries.  using Plots\npyplot() # Hide PyPlot warningsusing MIToS.MSA\n\nmsa = read(\"http://pfam.xfam.org/family/PF09776/alignment/full\", Stockholm)\n\nusing Plots\n\npyplot(size=(600,300))\n\nplot(   1:ncolumns(msa), # x is a range from 1 to the number of columns\n        vec(columngapfraction(msa)) .* 100.0, # y is a Vector{Float64} with the percentage of gaps of each column\n        linetype = :line,\n        ylabel = \"gaps [%]\",\n        xlabel = \"columns\",\n        legend=false)\n\npng(\"msa_gaps.png\") # hide\nnothing # hide(Image: )  plot(   1:nsequences(msa), # x is a range from 1 to the number of sequences\n        coverage(msa) .* 100, # y is a Vector{Float64} with the coverage of each sequence\n        linetype = :line,\n        ylabel = \"coverage [%]\",\n        xlabel = \"sequences\",\n        legend=false)\n\npng(\"msa_coverage.png\") # hide\nnothing # hide(Image: )  plot(msa)\npng(\"msa_msa.png\") # hide\nnothing # hide(Image: )  "
},

{
    "location": "MSA.html#Example:-Filter-sequences-per-coverage-and-columns-per-gap-fraction-1",
    "page": "MSA",
    "title": "Example: Filter sequences per coverage and columns per gap fraction",
    "category": "section",
    "text": "Taking advantage of the filter...! functions and the coverage and columngapfraction functions, it's possible to delete short sequences or columns with a lot of gaps.  println(\"\\tsequences\\tcolumns\")\nprintln( \"Before:\\t\", nsequences(msa), \"\\t\\t\", ncolumns(msa)  )\n# delete sequences with less than 90% coverage of the MSA length:\nfiltersequences!(msa, coverage(msa) .>= 0.9)\n# delete columns with more than 10% of gaps:\nfiltercolumns!(msa, columngapfraction(msa) .<= 0.1)\nprintln( \"After:\\t\", nsequences(msa), \"\\t\\t\",  ncolumns(msa)  )histogram(  vec(columngapfraction(msa)),\n            # Using vec() to get a Vector{Float64} with the fraction of gaps of each column\n            xlabel = \"gap fraction in [0,1]\", legend=false)\npng(\"msa_hist_gaps.png\") # hide\nnothing # hide(Image: )  histogram(  coverage(msa) .* 100.0, #  Column with the coverage of each sequence\n            xlabel = \"coverage [%]\", legend=false)\npng(\"msa_hist_coverage.png\") # hide\nnothing # hide(Image: )  "
},

{
    "location": "MSA.html#Example:-Plotting-the-percentage-of-identity-between-sequences-1",
    "page": "MSA",
    "title": "Example: Plotting the percentage of identity between sequences",
    "category": "section",
    "text": "The distribution of the percentage of identity between every pairs of sequences in a MSA, gives an idea of the MSA diversity. In this example, we are going to use percentidentity over a MSA to get that values.  using MIToS.MSA\nmsa = read(\"http://pfam.xfam.org/family/PF09776/alignment/full\", Stockholm)\npid = percentidentity(msa)\nnothing # hideMIToS stores the matrix of percentage of identity between the aligned sequences as a PairwiseListMatrix from the PairwiseListMatrices(Image: ) package. This matrix type saves RAM, allowing the storage of  big matrices. In this example, we use the to_table function of PairwiseListMatrices to convert the matrix into a table with indices.  using PairwiseListMatrices\n\npidtable = to_table(pid, diagonal=false)The function quantile gives a quick idea of the percentage identity distribution of the MSA.  quantile(convert(Vector{Float64}, pidtable[:,3]), [0.00, 0.25, 0.50, 0.75, 1.00])The function meanpercentidentity gives the mean value of the percent identity distribution for MSA with less than 300 sequences, or a quick estimate (mean PID in a random sample of sequence pairs) otherwise unless you set exact to true.  meanpercentidentity(msa)One can easily plot that matrix and its distribution using the heatmap and histogram functions of the Plots(Image: ) package.  using Plots\npyplot() # Hide PyPlot warningsusing Plots\npyplot()\nheatmap(full(pid), yflip=true, ratio=:equal)\npng(\"msa_heatmap_pid.png\") # hide\nnothing # hide(Image: )  histogram(pidtable[:,3], xlabel =\"Percentage of identity\", legend=false)\npng(\"msa_hist_pid.png\") # hide\nnothing # hide(Image: )  "
},

{
    "location": "MSA.html#Sequence-clustering-1",
    "page": "MSA",
    "title": "Sequence clustering",
    "category": "section",
    "text": "The MSA module allows to clusterize sequences in a MSA. The hobohmI function takes as input a MSA followed by an identity threshold value, and returns a Clusters type with the result of a Hobohm I(Image: ) sequence clustering. The Hobohm I algorithm will add a sequence to an existing cluster, if the percentage of identity is equal or greater than the threshold.   The Clusters is sub-type of ClusteringResult from the Clustering.jl(Image: ) package. One advantage of use a sub-type of ClusteringResultis that you are able to use any method defined on Clustering.jl like varinfo (Variation of Information) for example. Also, you can use any clustering algorithm included in Clustering.jl, and convert its result to an Clusters object to use them with MIToS.   MSA defines the functions nclusters to get the resulting number of clusters, counts to get the number of sequences on each cluster and assignments to get the cluster number of each sequence. The most important method is getweight, which returns the weight of each sequence. This method is used in the Information module of MIToS to reduce redundancy.  "
},

{
    "location": "MSA.html#Example:-Reducing-redundancy-of-a-MSA-1",
    "page": "MSA",
    "title": "Example: Reducing redundancy of a MSA",
    "category": "section",
    "text": "MSAs can suffer from an unnatural sequence redundancy and a high number of protein fragments. In this example, we will use a sequence clustering to make a non-redundant set of representative sequences using the function hobohmI to perform a clustering with the Hobohm I algorithm at 62% identity.  using Plots\npyplot() # Hide PyPlot warningsusing MIToS.MSA\n\nmsa = read(\"http://pfam.xfam.org/family/PF09776/alignment/full\", Stockholm)\n\nprintln(\"This MSA has \", nsequences(msa), \" sequences...\")clusters = hobohmI(msa, 62)println(\"...but has only \", nclusters(clusters), \" sequence clusters after a clustering at 62% identity.\")using Plots\npyplot()\n\nplot(msa)\npng(\"msa_clusters_i.png\") # hide\nnothing # hide(Image: )  We are going to use the DataFrames(Image: ) package to easily select the sequence with highest coverage of each cluster.  using DataFrames\n\ndf = DataFrame( seqnum = 1:nsequences(msa),\n                seqname = sequencenames(msa),\n                cluster = assignments(clusters), # the cluster number/index of each sequence\n                coverage = vec(coverage(msa)))It is possible to use this DataFrame and Plots to plot the sequence coverage of the MSA and also an histogram of the number of sequences in each cluster:  using StatPlots # Plotting DataFrames\nh = histogram(df[:cluster], ylabel=\"nseq\")\np = plot(df, :cluster, :coverage, linetype=:scatter)\nplot(p, h, nc=1, xlim=(0, nclusters(clusters)+1 ), legend=false)\npng(\"msa_clusters_ii.png\") # hide\nnothing # hide(Image: )  We use the Split-Apply-Combine strategy, though the by function of the DataFrames package, to select the sequence of highest coverage for each cluster.  maxcoverage = by(df, :cluster, cl -> cl[ findmax(cl[:coverage])[2] ,:])p = plot(maxcoverage, :cluster, :coverage, linetype=:scatter)\nh = histogram(maxcoverage[:cluster], ylabel=\"nseq\")\nplot(p, h, nc=1, xlim=(0, nclusters(clusters)+1 ), legend=false)\npng(\"msa_clusters_iii.png\") # hide\nnothing # hide(Image: )  We can easily generate a mask using list comprehension, to select only the representative sequences of the MSA (deleting the rest of the sequences with filtersequences!).  cluster_references = Bool[ seqnum in maxcoverage[:seqnum] for seqnum in 1:nsequences(msa) ]filtersequences!(msa, cluster_references)plot(msa)\npng(\"msa_clusters_iv.png\") # hide\nnothing # hide(Image: )  "
},

{
    "location": "Information.html#",
    "page": "Information",
    "title": "Information",
    "category": "page",
    "text": "CurrentModule = MIToS.Information"
},

{
    "location": "Information.html#Module-Information-1",
    "page": "Information",
    "title": "Information",
    "category": "section",
    "text": "The Information module of MIToS defines types and functions useful to calculate information measures (e.g. Mutual Information (MI) and Entropy) over a Multiple Sequence Alignment (MSA). This module was designed to count Residues (defined in the MSA module) in special contingency tables (as fast as possible) and to derive probabilities from these counts. Also, includes methods for applying corrections to those tables, e.g. pseudocounts and pseudo frequencies. Finally, Information allows to use these probabilities and counts to estimate information measures and other frequency based values.  using MIToS.Information # to load the Information module"
},

{
    "location": "Information.html#Features-1",
    "page": "Information",
    "title": "Features",
    "category": "section",
    "text": "Estimate multi dimensional frequencies and probability tables from sequences, MSAs, etc...\nCorrection for small number of observations\nCorrection for data redundancy on a MSA\nEstimate information measures\nCalculate corrected mutual information between residues  "
},

{
    "location": "Information.html#Contents-1",
    "page": "Information",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"Information.md\"]\nDepth = 4"
},

{
    "location": "Information.html#Counting-residues-1",
    "page": "Information",
    "title": "Counting residues",
    "category": "section",
    "text": "MIToS Information module defines a multidimensional ContingencyTable type and two types wrapping it, Counts and Probabilities, to store occurrences or probabilities. The ContingencyTable type stores the contingency matrix, its marginal values and total. These types are parametric, taking three ordered parameters:T : The type used for storing the counts or probabilities, e.g. Float64. It'spossible to use BigFloat if more precision it's needed.N : It's the dimension of the table and should be an Int.\nA : This should be a type, subtype of ResidueAlphabet, i.e.: UngappedAlphabet,GappedAlphabet or ReducedAlphabet.note: Note\nContingencyTable can be used for storing probabilities or counts. The wrapper types Probabilities and Counts are mainly intended to dispatch in methods that need to know if the matrix has probabilities or counts, e.g. entropy. In general, the use of ContingencyTable is recommended over the use of Probabilities and Counts.In this way, a matrix for storing pairwise probabilities of residues (without gaps) can be initialized using:  using MIToS.Information\n\nPij = ContingencyTable(Float64, Val{2}, UngappedAlphabet())[High level interface] It is possible to use the functions count and probabilities to easily calculate the frequencies of sequences or columns of a MSA, where the number of sequences/columns determine the dimension of the resulting table.  using MIToS.Information\nusing MIToS.MSA # to use res\"...\" to create Vector{Residue}\n\ncolumn_i = res\"AARANHDDRDC-\"\ncolumn_j = res\"-ARRNHADRAVY\"\n#   Nij[R,R] =   1     1   = 2\n\nNij = count(column_i, column_j)You can use sum to get the stored total:  sum(Nij) # There are 12 Residues, but 2 are gapsContingency tables can be indexed using Int or Residues:  Nij[2, 2] # Use Int to index the tableNij[Residue('R'), Residue('R')] # Use Residue to index the tablewarning: Warning\nThe number makes reference to the specific index in the table e.g [2,2] references the second row and the second column. The use of the number used to encode the residue to index the table is dangerous. The equivalent index number of a residue depends on the used alphabet and Int(Residue('X')) will be always out of bounds.  Indexing with Residues works as expected. It uses the alphabet of the contingency table to find the index of the Residue.using MIToS.Information\nusing MIToS.MSA\n\nalphabet = ReducedAlphabet(\"(AILMV)(NQST)(RHK)(DE)(FWY)CGP\")\n\ncolumn_i = res\"AARANHDDRDC-\"\ncolumn_j = res\"-ARRNHADRAVY\"\n#   Fij[R,R] =   1  1  1   = 3 # RHK\n\nFij = count(column_i, column_j, alphabet=alphabet)Fij[Residue('R'), Residue('R')] # Use Residue to index the tableThe function getcontingencytable allows to access the wrapped ContingencyTable in a Counts object. You can use it, in combination with normalize to get a contingency table of probabilities. The result can be wrapped inside a Probabilities object:  Probabilities(normalize(getcontingencytable(Fij)))"
},

{
    "location": "Information.html#Example:-Plotting-the-probabilities-of-each-residue-in-a-sequence-1",
    "page": "Information",
    "title": "Example: Plotting the probabilities of each residue in a sequence",
    "category": "section",
    "text": "Similar to the count function, the probabilities function can take at least one sequence (vector of residues) and returns the probabilities of each residue. Optionally, the keyword argument alphabet could be used to count some residues in the same cell of the table.  probabilities(res\"AARANHDDRDC\", alphabet=alphabet)Here, we are going to use the probabilities function to get the residue probabilities of a particular sequence from UniProt.use the getsequence function, from the MSA module, to get the sequence from a FASTA downloaded from UniProt.  using MIToS.Information # to use the probabilities function\nusing MIToS.MSA # to use getsequence on the one sequence FASTA (canonical) from UniProt\nseq = read(\"http://www.uniprot.org/uniprot/P29374.fasta\", FASTA) # Small hack: read the single sequence as a MSA\nprobabilities(seq[1,:]) # Select the single sequence and calculate the probabilitiesnote: Note\nIn the previous example, using getsequence(seq,1) instead of seq[1,:] will return the sequence as a matrix with a single column to keep information for both dimensions. To use probabilities (or count) you can make use of the Julia's vec function to transform the matrix to a vector, e.g.: probabilities(vec(getsequence(seq,1))).using Plots\npyplot(size=(600,300))\nusing MIToS.Information # to use the probabilities function\nusing MIToS.MSA # to use getsequence on the one sequence FASTA (canonical) from UniProt\nseq = read(\"http://www.uniprot.org/uniprot/P29374.fasta\", FASTA) # Small hack: read the single sequence as a MSA\nfrequencies = probabilities(seq[1,:]) # Select the single sequence and calculate the probabilitiesusing Plots # We choose Plots because it's intuitive, concise and backend independent\npyplot(size=(600,300))You can plot together with the probabilities of each residue in a given sequence, the probabilities of each residue estimated with the BLOSUM62 substitution matrix. That matrix is exported as a constant by the Information module as BLOSUM62_Pi.  bar(\n    1:20,\n    [ frequencies  BLOSUM62_Pi ],\n    lab = [ \"Sequence\"  \"BLOSUM62\"   ],\n    alpha=0.5\n    )\npng(\"inf_plotfreq.png\") # hide\nnothing # hide(Image: )  "
},

{
    "location": "Information.html#Low-count-corrections-1",
    "page": "Information",
    "title": "Low count corrections",
    "category": "section",
    "text": "Low number of observations can lead to sparse contingency tables, that lead to wrong probability estimations. It is shown in Buslje et. al. 2009(Image: ) that low-count corrections, can lead to improvements in the contact prediction capabilities of the Mutual Information. The Information module has available two low-count corrections:  Additive Smoothing(Image: ); the constant value pseudocount described in Buslje et. al. 2009(Image: ).  \nBLOSUM62 based pseudo frequencies of residues pairs, similar to Altschul et. al. 1997(Image: ).  using MIToS.MSA\n\nmsa = read(\"http://pfam.xfam.org/family/PF09776/alignment/full\", Stockholm)\n\nfiltercolumns!(msa, columngapfraction(msa) .< 0.5) # delete columns with 50% gaps or more\n\ncolumn_i = msa[:,1]\ncolumn_j = msa[:,2]If you have a preallocated ContingencyTable you can use count! to fill it, this prevent to create a new table as count do. However, you should note that count! adds the new counts to the pre existing values, so in this case, we want to start with a table initialized with zeros.  using MIToS.Information\n\nconst alphabet = ReducedAlphabet(\"(AILMV)(NQST)(RHK)(DE)(FWY)CGP\")\n\nNij = ContingencyTable(Float64, Val{2}, alphabet)#      table  weights         pseudocount      sequences...\ncount!(Nij,   NoClustering(), NoPseudocount(), column_i, column_j)note: Note\nYou can use NoClustering() in places where clustering weights are required to not use weights. Also, NoPseudocount() in places where pseudocount values are required to not use pseudocounts.In cases like the above, where there are few observations, it is possible to apply a constant pseudocount to the counting table.  This module defines the type AdditiveSmoothing and the correspond fill! and  apply_pseudocount! methods to efficiently add or fill with a constant value each element of the table.apply_pseudocount!(Nij, AdditiveSmoothing(1.0))[High level interface.] The count function has a pseudocounts keyword argument that can take a AdditiveSmoothing value to easily calculate occurrences with pseudocounts. Also the alphabet keyword argument can be used to chage the default alphabet (i.e. )count(column_i, column_j, pseudocounts=AdditiveSmoothing(1.0), alphabet=alphabet)To use the conditional probability matrix BLOSUM62_Pij in the calculation of pseudo frequencies G for the pair of residues a, b, it should be calculated first the real frequencies/probabilities p_ab. The observed probabilities are then used to estimate the pseudo frequencies.  G_ab = sum_cd  p_cd cdot BLOSUM62( a  c ) cdot BLOSUM62( b  d )Finally, the probability P of each pair of residues a, b between the columns i, j is the weighted mean between the observed frequency p and BLOSUM62-based pseudo frequency G, where α is generally the number of clusters or the number of sequences of the MSA and β is an empiric weight value. β was determined to be close to 8.512.  P_ab = fracalpha cdot p_ab + beta cdot G_ab alpha + betaThis could be easily achieved using the pseudofrequencies keyword argument of the probabilities function. That argument can take a BLOSUM_Pseudofrequencies object that is created with α and β as first and second argument, respectively.Pij = probabilities(column_i, column_j, pseudofrequencies=BLOSUM_Pseudofrequencies(nsequences(msa), 8.512))You can also use apply_pseudofrequencies! in a previously filled probability contingency table. i.e. apply_pseudofrequencies!(Pij, BLOSUM_Pseudofrequencies(α, β))warning: Warning\nBLOSUM_Pseudofrequencies can be only be applied in normalized/probability tables with UngappedAlphabet.  "
},

{
    "location": "Information.html#Correction-for-data-redundancy-in-a-MSA-1",
    "page": "Information",
    "title": "Correction for data redundancy in a MSA",
    "category": "section",
    "text": "A simple way to reduce redundancy in a MSA without losing sequences, is clusterization and sequence weighting. The weight of each sequence should be 1/N, where N is the number of sequences in its cluster. The Clusters type of the MSA module stores the weights. This vector of weights can be extracted (with the getweight function) and used by the count and probabilities functions with the keyword argument weights. Also it's possible to use the Clusters as second argument of the function count!.  clusters = hobohmI(msa, 62) # from MIToS.MSAcount(msa[:,1], msa[:,2], weights=clusters)"
},

{
    "location": "Information.html#Estimating-information-measures-on-an-MSA-1",
    "page": "Information",
    "title": "Estimating information measures on an MSA",
    "category": "section",
    "text": "The Information module has a number of functions defined to calculate information measures from Counts and Probabilities:entropy : Shannon entropy (H)\nmarginal_entropy : Shannon entropy (H) of the marginals\nkullback_leibler : Kullback-Leibler (KL) divergence\nmutual_information : Mutual Information (MI)\nnormalized_mutual_information : Normalized Mutual Information (nMI) by Entropy\ngap_intersection_percentage\ngap_union_percentageInformation measure functions take optionally the base as the last positional argument (default: e). You can use 2.0 to measure information in bits.using MIToS.Information\nusing MIToS.MSA\n\nNi = count(res\"PPCDPPPPPKDKKKKDDGPP\") # Ni has the count table of residues in this low complexity sequence\n\nH = entropy(Ni) # returns the Shannon entropy in nats (base e)H = entropy(Ni, 2.0) # returns the Shannon entropy in bits (base 2)Information module defines special iteration functions to easily and efficiently compute a measure over a MSA. In particular, mapcolfreq! and mapseqfreq! map a function that takes a table of Counts or Probabilities. The table is filled in place with the counts or probabilities of each column or sequence of a MSA, respectively. mapcolpairfreq! and mapseqpairfreq! are similar, but they fill the table using pairs of columns or sequences, respectively.  This functions take three positional arguments: the function f to be calculated, the msa and table of Counts or Probabilities.  After that, this function takes some keyword arguments:weights (default: NoClustering()) : Weights to be used for table counting.\npseudocounts (default: NoPseudocount()) : Pseudocount object to be applied to table.\npseudofrequencies (default: NoPseudofrequencies()) : Pseudofrequencies to beapplied to the normalized (probabilities) table.  mapcolpairfreq! and mapseqpairfreq! also have a fourth positional argument usediagonal that indicates if the function should be applied to identical element pairs (default to Val{true}). This two functions also have an extra keyword argument diagonalvalue (default to zero) to indicate the value used to fill the diagonal elements if usediagonal is Val{false}.  "
},

{
    "location": "Information.html#Example:-Estimating-*H(X)*-and-*H(X,-Y)*-over-an-MSA-1",
    "page": "Information",
    "title": "Example: Estimating H(X) and H(X, Y) over an MSA",
    "category": "section",
    "text": "In this example, we are going to use mapcolfreq! and mapcolpairfreq! to estimate Shannon entropy of MSA columns H(X) and the joint entropy H(X, Y) of columns pairs, respectively.  using Plots\npyplot()using MIToS.MSA\n\nmsa = read(\"http://pfam.xfam.org/family/PF09776/alignment/full\", Stockholm)We are going to count residues to estimate the entropy. The entropy estimation is performed over a rehused Counts object. The result will be a vector containing the values estimated over each column without counting gaps (UngappedAlphabet).  using MIToS.Information\n\nHx = mapcolfreq!(entropy, msa, Counts(ContingencyTable(Float64, Val{1}, UngappedAlphabet())))If we want the joint entropy between columns pairs, we need to use a bidimensional table of Counts and mapcolpairfreq!.Hxy = mapcolpairfreq!(entropy, msa, Counts(ContingencyTable(Float64, Val{2}, UngappedAlphabet())))In the above examples, we indicate the type of each occurrence in the counting and the probability table to use. Also, it's possible for some measures as entropy and mutual information, to estimate the values only with the count table (without calculate the probability table). Estimating measures only with a ResidueCount table, when this is possible, should be faster than using a probability table.  Time_Pab = map(1:100) do x\n    time = @elapsed mapcolpairfreq!(entropy, msa, Probabilities(ContingencyTable(Float64, Val{2}, UngappedAlphabet())))\nend\n\nTime_Nab = map(1:100) do x\n    time = @elapsed mapcolpairfreq!(entropy, msa, Counts(ContingencyTable(Float64, Val{2}, UngappedAlphabet())))\nend\n\nusing Plots\npyplot()\n\nhistogram( [Time_Pab Time_Nab],\n    labels = [\"Using ResidueProbability\" \"Using ResidueCount\"],\n    xlabel = \"Execution time [seconds]\" )\n\npng(\"inf_entropy.png\") # hide\nnothing # hide(Image: )   "
},

{
    "location": "Information.html#MIToS.Information.buslje09",
    "page": "Information",
    "title": "MIToS.Information.buslje09",
    "category": "Function",
    "text": "buslje09 takes a MSA or a file and a Format as first arguments. It calculates a Z score and a corrected MI/MIp as described on Busjle et. al. 2009.\n\nkeyword argument, type, default value and descriptions:\n\n  - lambda      Float64   0.05    Low count value\n  - clustering  Bool      true    Sequence clustering (Hobohm I)\n  - threshold             62      Percent identity threshold for clustering\n  - maxgap      Float64   0.5     Maximum fraction of gaps in positions included in calculation\n  - apc         Bool      true    Use APC correction (MIp)\n  - samples     Int       100     Number of samples for Z-score\n  - fixedgaps   Bool      true    Fix gaps positions for the random samples\n  - alphabet    ResidueAlphabet UngappedAlphabet()  Residue alphabet to be used\n\nThis function returns:\n\n  - Z score\n  - MI or MIp\n\n\n\n"
},

{
    "location": "Information.html#MIToS.Information.BLMI",
    "page": "Information",
    "title": "MIToS.Information.BLMI",
    "category": "Function",
    "text": "BLMI takes a MSA or a file and a Format as first arguments. It calculates a Z score (ZBLMI) and a corrected MI/MIp as described on Busjle et. al. 2009 but using using BLOSUM62 pseudo frequencies instead of a fixed pseudocount.\n\nKeyword argument, type, default value and descriptions:\n\n  - beta        Float64   8.512   β for BLOSUM62 pseudo frequencies\n  - lambda      Float64   0.0     Low count value\n  - threshold             62      Percent identity threshold for sequence clustering (Hobohm I)\n  - maxgap      Float64   0.5     Maximum fraction of gaps in positions included in calculation\n  - apc         Bool      true    Use APC correction (MIp)\n  - samples     Int       50      Number of samples for Z-score\n  - fixedgaps   Bool      true    Fix gaps positions for the random samples\n\nThis function returns:\n\n  - Z score (ZBLMI)\n  - MI or MIp using BLOSUM62 pseudo frequencies (BLMI/BLMIp)\n\n\n\n"
},

{
    "location": "Information.html#Corrected-Mutual-Information-1",
    "page": "Information",
    "title": "Corrected Mutual Information",
    "category": "section",
    "text": "MIToS ships with two methods to easily calculate corrected mutual information.   The first is the algorithm described in Buslje et. al. 2009(Image: ). This algorithm can be accessed through the buslje09 function and includes:  Low count correction using AdditiveSmoothing\nSequence weighting after a hobohmI clustering\nAverage Product Correction (APC) proposed byDunn et. al. 2008(Image: ), through the APC! function that takes a MI matrix.Z score correction using the functions shuffle! from the MSA module and zscorefrom the PairwiseListMatrices package.  buslje09The second, implemented in the BLMI function, has the same corrections that the above algorithm, but use BLOSUM62 pseudo frequencies. This function is slower than buslje09 (at the same number of samples), but gives better performance (for structural contact prediction) when the MSA has less than 400 clusters after a Hobohm I at 62% identity.  BLMI"
},

{
    "location": "Information.html#Example:-Estimating-corrected-MI-from-an-MSA-1",
    "page": "Information",
    "title": "Example: Estimating corrected MI from an MSA",
    "category": "section",
    "text": "using Plots\npyplot()using MIToS.MSA\nusing MIToS.Information\n\nmsa = read(\"http://pfam.xfam.org/family/PF16078/alignment/full\", Stockholm)\nZMIp, MIp  = buslje09(msa)\nZMIpZBLMIp, BLMIp  = BLMI(msa)\nZBLMIp"
},

{
    "location": "Information.html#Visualize-Mutual-Information-1",
    "page": "Information",
    "title": "Visualize Mutual Information",
    "category": "section",
    "text": "You can use the function of the Plots package to visualize the Mutual Information (MI) network between residues. As an example, we are going to visualize the MI between residues of the Pfam domain PF16078. The heatmap is the simplest way to visualize the values of the Mutual Information matrix.  using Plots\npyplot()\n\nheatmap(ZMIp, yflip=true)\npng(\"inf_heatmap.png\") # hide\nnothing # hide(Image: )   ZMIp is a Z score of the corrected MIp against its distribution on a random MSA (shuffling the residues in each sequence), so pairs with highest values are more likely to co-evolve. Here, we are going to use the top 1% pairs of MSA columns.  using PairwiseListMatrices # to use getlist\n\nthreshold = quantile(getlist(ZMIp), 0.99)ZMIp[ ZMIp .< threshold ] = NaN\nheatmap(ZMIp, yflip=true)\npng(\"inf_heatmap_top.png\") # hide\nnothing # hide(Image: )   We are going to calculate the cMI (cumulative mutual information) value of each node. Where cMI is a mutual information score per position that characterizes the extent of mutual information \"interactions\" in its neighbourhood. This score is calculated as the sum of MI values above a certain threshold for every amino acid pair where the particular residue appears. This value defines to what degree a given amino acid takes part in a mutual information network and we are going to indicate it using the node color. To calculate cMI we are going to use the cumulative function:   cMI = cumulative(ZMIp, threshold)The nodes have an order, because they are columns in a MSA. So, the arc diagram it's useful to visualize long and short association between MSA positions. In general, long interactions has more interest.using PlotRecipes\n\ngraphplot(ZMIp, size=(600,250), method=:arcdiagram) # , zcolor=cMI)You can also use a chord diagram to see the same pattern.  graphplot(ZMIp, size=(600,600), method=:chorddiagram)"
},

{
    "location": "SIFTS.html#",
    "page": "SIFTS",
    "title": "SIFTS",
    "category": "page",
    "text": ""
},

{
    "location": "SIFTS.html#Module-SIFTS-1",
    "page": "SIFTS",
    "title": "SIFTS",
    "category": "section",
    "text": "The SIFTS module of MIToS allows to obtain the residue-level mapping between databases stored in the SIFTS XML files. It makes easy to assign PDB residues to UniProt/Pfam positions.   Given the fact that pairwise alignments can lead to misleading association between residues in both sequences, SIFTS offers  more reliable association between sequence and structure residue numbers.  using MIToS.SIFTS # to load the SIFTS module"
},

{
    "location": "SIFTS.html#Features-1",
    "page": "SIFTS",
    "title": "Features",
    "category": "section",
    "text": "Download and parse SIFTS XML files\nStore residue-level mapping in Julia\nEasy generation of Dicts between residues numbers"
},

{
    "location": "SIFTS.html#Contents-1",
    "page": "SIFTS",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"SIFTS.md\"]\nDepth = 4"
},

{
    "location": "SIFTS.html#Simplest-residue-level-mapping-1",
    "page": "SIFTS",
    "title": "Simplest residue-level mapping",
    "category": "section",
    "text": "This module export the function siftsmapping to generate a Dict between residue numbers. This function takes 5 positional arguments.     1) The name of the SIFTS XML file to parse,       2) the source database,       3) the source protein/structure identifier,       4) the destiny database and,       5) the destiny protein/structure identifier.   Optionally it’s possible to indicate a particular PDB chain and if missings will be used.  Databases should be indicated using an available sub-type of DataBase. Keys and values types will be depend on the residue number type in that database.Type db... Database Residue number type\ndbPDBe PDBe (Protein Data Bank in Europe) Int\ndbInterPro InterPro ASCIIString\ndbUniProt UniProt Int\ndbPfam Pfam (Protein families database) Int\ndbNCBI NCBI (National Center for Biotechnology Information) Int\ndbPDB PDB (Protein Data Bank) ASCIIString\ndbCATH CATH ASCIIString\ndbSCOP SCOP (Structural Classification of Proteins) ASCIIStringTo download the XML SIFTS file of a determined PDB use the downloadsifts function.  using MIToS.SIFTS\n\nsiftsfile = downloadsifts(\"1IVO\")The following example, shows the residue number mapping between Pfam and PDB. Pfam uses UniProt coordinates and PDB uses their own residue numbers with insertion codes. Note that the siftsmapping function is case sensitive, and that SIFTS stores PDB identifiers using lowercase characters.  siftsmap = siftsmapping(siftsfile,\n                        dbPfam,\n                        \"PF00757\",\n                        dbPDB,\n                        \"1ivo\", # SIFTS stores PDB identifiers in lowercase\n                        chain=\"A\", # In this example we are only using the chain A of the PDB\n                        missings=false) # Residues without coordinates aren't used in the mapping"
},

{
    "location": "SIFTS.html#Storing-residue-level-mapping-1",
    "page": "SIFTS",
    "title": "Storing residue-level mapping",
    "category": "section",
    "text": "If you need more than the residue number mapping between two databases, you could access all the residue-level cross references using the function read in the SIFTSXMLFormat file. The parse function (and therefore the read function) for the SIFTSXML format, also takes the keyword arguments chain and missings. The read/parse function returns a Vector of SIFTSResidues objects that stores the cross references between residues in each database.  siftsresidues = read(siftsfile, SIFTSXML, chain=\"A\", missings=false) # Array{SIFTSResidue,1}\n\nresidue_data = siftsresidues[300]You are free to access the SIFTSResidue fields in order to get the desired information. SIFTSResidue objects contain db... objects (sub-types of DataBase), with the cross referenced information. You should note that, except for the PDBe and InterPro fields, the fields are Nullables objects so, you need to use the get function to access the db... object. For example, getting the UniProt residue name (one letter code of the amino acid) would be:  isnull(residue_data.UniProt) ? \"\" : get(residue_data.UniProt).nameThat line of code returns an empty string if the UniProt field is null. Otherwise, it returns a string with the name of the residue in UniProt. Because that way of access values in a Residue is too verbose, MIToS defines a more complex signature for get. Using MIToS get the previous line of code will be:  #   SIFTSResidue  database   field  default\nget(residue_data, dbUniProt, :name, \"\")The is not need to use the full signature, but the returned value will change. In particular, a Nullable object is returned if a default value is not given at the end of the signature:  using MIToS.SIFTS; residue_data = read(\"1ivo.xml.gz\", SIFTSXML)[301] # hide\nget(residue_data, dbUniProt) # Takes the database type and returns a nullable with the field content\nget(residue_data, dbUniProt, :name) # Takes also a Symbol with a field name and returns a nullable with the field content inside the database typeBut you don't need the getfunction to access the three letter code of the residue in PDBe because the PDBe field is not Nullable.residue_data.PDBe.nameSIFTSResidue also store information about if that residue is missing in the PDB structure and the information about the secondary structure (sscode and ssname):  using MIToS.SIFTS; residue_data = read(\"1ivo.xml.gz\", SIFTSXML)[301] # hide\nresidue_data.missing\nresidue_data.sscode\nresidue_data.ssname"
},

{
    "location": "SIFTS.html#Accessing-residue-level-cross-references-1",
    "page": "SIFTS",
    "title": "Accessing residue-level cross references",
    "category": "section",
    "text": "You can ask for particular values in a single SIFTSResidue using the get function.  using MIToS.SIFTS\nresidue_data = read(\"1ivo.xml.gz\", SIFTSXML)[301]\n# Is the UniProt residue name in the list of basic amino acids [\"H\", \"K\", \"R\"]?\nget(residue_data, dbUniProt, :name, \"\") in [\"H\", \"K\", \"R\"]Use higher order functions and lambda expressions (anonymous functions) or list comprehension to easily ask for information on the Vector{SIFTSResidue}. You can use get with the previous signature or simple get, direct field access and isnull.# Captures PDB residue numbers if the Pfam id is \"PF00757\"\nresnums = [ get(res.PDB).number for res in siftsresidues if !isnull(res.PDB) && get(res, dbPfam, :id, \"\") == \"PF00757\" ]Useful higher order functions are:find  # Which of the residues have UniProt residue names in the list [\"H\", \"K\", \"R\"]? (basic residues)\nindexes = find(res -> get(res, dbUniProt, :name, \"\") in [\"H\", \"K\", \"R\"], siftsresidues)map  map(i -> get(siftsresidues[i].UniProt), indexes) # UniProt data of the basic residuesfilter  # SIFTSResidues with UniProt names in [\"H\", \"K\", \"R\"]\nbasicresidues = filter(res -> get(res, dbUniProt, :name, \"\") in [\"H\", \"K\", \"R\"], siftsresidues)\n\nget(basicresidues[1].UniProt) # UniProt data of the first basic residue"
},

{
    "location": "SIFTS.html#Example:-Which-residues-are-missing-in-the-PDB-structure-1",
    "page": "SIFTS",
    "title": "Example: Which residues are missing in the PDB structure",
    "category": "section",
    "text": "Given that SIFTSResidue objects store a missing residue flag, it’s easy to get a vector where there is a true value if the residue is missing in the structure.  using MIToS.SIFTS\nsifts_1ivo = read(\"1ivo.xml.gz\", SIFTSXML, chain=\"A\"); # SIFTSResidues of the 1IVO chain A\n[res.missing for res in sifts_1ivo]However, if you need to filter using other conditions, you’ll find useful the get function. In this example, we are going to ask for the UniProt id (to avoid problems with fragments, tags or chimeric/fusion proteins). We are also using get to select an specific PDB chain.  using MIToS.SIFTS\ndownloadsifts(\"1JQZ\")using MIToS.SIFTS\nsifts_1jqz = read(\"1jqz.xml.gz\", SIFTSXML); # It has an amino terminal his tag\nmissings = [ (  ( get(res, dbUniProt, :id, \"\") == \"P05230\" ) &\n                ( get(res, dbPDB, :chain, \"\") ==  \"A\" ) &\n                res.missing ) for res in sifts_1jqz             ];\nprintln(\"There are only \", sum(missings), \" missing residues in the chain A, associated to UniProt P05230\")\nprintln(\"But there are \", sum([ res.missing for res in sifts_1jqz ]), \" missing residues in the PDB file.\")"
},

{
    "location": "PDB.html#",
    "page": "PDB",
    "title": "PDB",
    "category": "page",
    "text": ""
},

{
    "location": "PDB.html#Module-PDB-1",
    "page": "PDB",
    "title": "PDB",
    "category": "section",
    "text": "The module PDB defines types and methods to work with protein structures inside Julia. It is useful to link structural and sequential information, and needed for measure the predictive performance at protein contact prediction of mutual information scores.  using MIToS.PDB # to load the PDB module"
},

{
    "location": "PDB.html#Features-1",
    "page": "PDB",
    "title": "Features",
    "category": "section",
    "text": "Read and parse PDB and PDBML files.\nCalculate distance and contacts between atoms or residues.\nDetermine interaction between residues."
},

{
    "location": "PDB.html#Contents-1",
    "page": "PDB",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"PDB.md\"]\nDepth = 4"
},

{
    "location": "PDB.html#Retrieve-information-from-PDB-database-1",
    "page": "PDB",
    "title": "Retrieve information from PDB database",
    "category": "section",
    "text": "This module exports the downloadpdb function, to retrieve a PDB file from   PDB database(Image: ). This function downloads a gzipped PDBML file, which could be easily read it with MIToS by default, but you are able to determine the format as PDBFile if you want it.  using MIToS.PDB\n\npdbfile = downloadpdb(\"1IVO\", format=PDBFile)PDB module also exports a getpdbdescription to access the header information of a PDB entry.  getpdbdescription(\"1IVO\")"
},

{
    "location": "PDB.html#Read-and-parse-PDB-files-1",
    "page": "PDB",
    "title": "Read and parse PDB files",
    "category": "section",
    "text": "This is easy using the read and parse functions, indicating the filename and the Format: PDBML for PDB XML files or PDBFile for usual PDB files. These functions returns a Vector of PDBResidue objects with all the residues in the PDB.   To return only a specific subset of residues/atoms you can use any of the following keyword arguments:  keyword arguments default returns only ...\nchain All residues from a PDB chain, i.e. \"A\"\nmodel All residues from a determined model, i.e. \"1\"\ngroup All residues from a group: \"ATOM\", \"HETATM\" or All for both\natomname All atoms with a specific name, i.e. \"CA\"\nonlyheavy false heavy atoms (not hydrogens) if it's true\noccupancyfilter false only the atoms with the best occupancy are returned if it's truenote: Note\nFor PDBML files it is possible to use the keyword argument label to false (default to true) to get the auth_ attributes instead of the label_ attributes for chain, atom and residue name fields. The auth_ attributes are alternatives provided by an author in order to match the identification/values used in the publication that describes the structure.  # Read α carbon of each residue from the 1ivo pdb file, in the model 1, chain A and in the ATOM group.\nCA_1ivo = read(pdbfile, PDBFile, model=\"1\", chain=\"A\", group=\"ATOM\", atomname=\"CA\")\n\nCA_1ivo[1] # First residue. It has only the α carbon."
},

{
    "location": "PDB.html#Looking-for-particular-residues-1",
    "page": "PDB",
    "title": "Looking for particular residues",
    "category": "section",
    "text": "MIToS parse PDB files to vector of residues, instead of using a hierarchical structure like other packages. This approach makes the search and selection of residues or atoms a little different. To make it easy, this module exports a number of functions and macros to select particular residues or atoms. Given the fact that residue numbers from different chains, models, etc. can collide, it's mandatory to indicate the model, chain, group, residue number and atom name in a explicit way to these functions or macros. If you want to select all the residues in one of the categories, you are able to use the type All. You can also use regular expressions or functions to make the selections.using MIToS.PDB\npdbfile = downloadpdb(\"1IVO\", format=PDBFile)\nresidues_1ivo = read(pdbfile, PDBFile)\n# Select residue number 9 from model 1 and chain B\nresidues(residues_1ivo, \"1\", \"B\", All, \"9\")"
},

{
    "location": "PDB.html#Getting-a-Dict-of-PDBResidues-1",
    "page": "PDB",
    "title": "Getting a Dict of PDBResidues",
    "category": "section",
    "text": "If you prefer a Dict of PDBResidue, indexed by their residue numbers, you can use the residuedict function or the @residuedict macro.  # Dict of residues from the model 1, chain A and from the ATOM group\nchain_a = residuesdict(residues_1ivo, \"1\", \"A\", \"ATOM\", All)\nchain_a[\"9\"]You can do the same with the macro @residuesdict to get a more readable code  chain_a = @residuesdict residues_1ivo model \"1\" chain \"A\" group \"ATOM\" residue All\nchain_a[\"9\"]"
},

{
    "location": "PDB.html#Select-particular-residues-1",
    "page": "PDB",
    "title": "Select particular residues",
    "category": "section",
    "text": "Use the residues function to collect specific residues. It's possible to use a single residue number (i.e. \"2\") or even a function which should return true for the selected residue numbers. Also regular expressions can be used to select residues. Use All to select all the residues.  residue_list = map(string, 2:5)\n\n# If the list is large, you can use a `Set` to gain performance\n# residue_set = Set(map(string, 2:5))first_res = residues(residues_1ivo, \"1\", \"A\", \"ATOM\", resnum -> resnum in residue_list)\n\nfor res in first_res\n    println(res.id.name, \" \", res.id.number)\nendA more complex example using an anonymous function:  # Select all the residues of the model 1, chain A of the ATOM group with residue number less than 5\n\nfirst_res = residues(residues_1ivo, \"1\", \"A\", \"ATOM\", x -> parse(Int, match(r\"^(\\d+)\", x)[1]) <= 5 )\n# The anonymous function takes the residue number (string) and use a regular expression\n# to extract the number (without insertion code).\n# It converts the number to `Int` to test if the it is `<= 5`.\n\nfor res in first_res\n    println(res.id.name, \" \", res.id.number)\nendUse the @residues macro for a cleaner syntax.  # You can use All, regular expressions or functions also for model, chain and group:\n\n# i.e. Takes the residue 10 from chains A and B\n\nfor res in @residues residues_1ivo model \"1\" chain ch -> ch in [\"A\",\"B\"] group \"ATOM\" residue \"10\"\n    println(res.id.chain, \" \", res.id.name, \" \", res.id.number)\nend"
},

{
    "location": "PDB.html#Select-particular-atoms-1",
    "page": "PDB",
    "title": "Select particular atoms",
    "category": "section",
    "text": "The atoms function or macro allow to select a particular set of atoms.# Select all the atoms with name starting with \"C\" using a regular expression\n# from all the residues of the model 1, chain A of the ATOM group\n\ncarbons = @atoms residues_1ivo model \"1\" chain \"A\" group \"ATOM\" residue All atom r\"C.+\"\n\ncarbons[1]You can also use the atoms function instead of the @atoms macro:  atoms(residues_1ivo, \"1\", \"A\", \"ATOM\", All, r\"C.+\")[1]"
},

{
    "location": "PDB.html#Protein-contact-map-1",
    "page": "PDB",
    "title": "Protein contact map",
    "category": "section",
    "text": "The PDB module offers a number of functions to measure distances between atoms or residues, to detect possible interactions or contacts. In particular the contact function calls the distance function using a threshold or limit in an optimized way. The measure can be done between alpha carbons (\"CA\"), beta carbons (\"CB\") (alpha carbon for glycine), any heavy atom (\"Heavy\") or any (\"All\") atom of the residues.In the following example, whe are going to plot a contact map for the 1ivo chain A. Two residues will be considered in contact if their β carbons (α carbon for glycine) have a distance of 8Å or less.  using MIToS.PDB\n\npdbfile = downloadpdb(\"1IVO\", format=PDBFile)\n\nresidues_1ivo = read(pdbfile, PDBFile)\n\npdb = @residues residues_1ivo model \"1\" chain \"A\" group \"ATOM\" residue All\n\ndmap = distance(pdb, criteria=\"All\") # Minimum distance between residues using all their atomsUse the contact function to get a contact map:  cmap = contact(pdb, 8.0, criteria=\"CB\") # Contact mapusing Plots\npyplot() # Hide PyPlot warningsusing Plots\npyplot()\n\nheatmap(dmap, grid=false, yflip=true, ratio=:equal)\n\npng(\"pdb_dmap.png\") # hide\nnothing # hide(Image: )  heatmap(cmap, grid=false, yflip=true, ratio=:equal)\n\npng(\"pdb_cmap.png\") # hide\nnothing # hide(Image: )  "
},

{
    "location": "PDB.html#Structural-superposition-1",
    "page": "PDB",
    "title": "Structural superposition",
    "category": "section",
    "text": "using Plots\npyplot() # Hide PyPlot warningsusing MIToS.PDB\n\npdbfile = downloadpdb(\"2HHB\")\n\nres_2hhb = read(pdbfile, PDBML)\n\nchain_A = pdb = @residues res_2hhb model \"1\" chain \"A\" group \"ATOM\" residue All\nchain_C = pdb = @residues res_2hhb model \"1\" chain \"C\" group \"ATOM\" residue All\n\nusing Plots\npyplot()\n\nscatter3d(chain_A, label=\"A\", alpha=0.5)\nscatter3d!(chain_C, label=\"C\", alpha=0.5)\n\npng(\"pdb_unaligned.png\") # hide\nnothing # hide(Image: )  superimposed_A, superimposed_C, RMSD = superimpose(chain_A, chain_C)\n\nRMSDscatter3d(superimposed_A, label=\"A\", alpha=0.5)\nscatter3d!(superimposed_C, label=\"C\", alpha=0.5)\npng(\"pdb_aligned.png\") # hide\nnothing # hide(Image: )  "
},

{
    "location": "Pfam.html#",
    "page": "Pfam",
    "title": "Pfam",
    "category": "page",
    "text": ""
},

{
    "location": "Pfam.html#Module-Pfam-1",
    "page": "Pfam",
    "title": "Pfam",
    "category": "section",
    "text": "MIToS defines methods and types useful for any MSA. The Pfam module uses other MIToS modules in the context of Pfam MSAs, where it’s possible to us determine how structure and sequence information should be mapped. This module defines functions that go from a Pfam MSA to the protein contact prediction performance of pairwise scores estimated from that MSA.using MIToS.Pfam # to load the Pfam module"
},

{
    "location": "Pfam.html#Features-1",
    "page": "Pfam",
    "title": "Features",
    "category": "section",
    "text": "Download and read Pfam MSAs.\nObtain PDB information from alignment annotations.\nMap between sequence/alignment residues/columns and PDB structures.\nMeasure of AUC (ROC curve) for protein contact prediction of MI scores."
},

{
    "location": "Pfam.html#Contents-1",
    "page": "Pfam",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"Pfam.md\"]\nDepth = 4"
},

{
    "location": "Pfam.html#Getting-a-Pfam-MSA-1",
    "page": "Pfam",
    "title": "Getting a Pfam MSA",
    "category": "section",
    "text": "The function downloadpfam takes a Pfam accession and downloads a Pfam MSA in Stockholm format. Use read function and the Stockholm Format to get a AnnotatedMultipleSequenceAlignment object with the MSA and its Pfam annotations. You must set generatemapping and useidcoordinates to true the first time you read the downloaded MSA. This is necessary to some of the methods in the Pfam module.  using MIToS.Pfam\npfamfile = downloadpfam(\"PF12464\")\nmsa = read(pfamfile, Stockholm, generatemapping=true, useidcoordinates=true)"
},

{
    "location": "Pfam.html#Getting-PDB-information-from-an-MSA-1",
    "page": "Pfam",
    "title": "Getting PDB information from an MSA",
    "category": "section",
    "text": "The function getseq2pdb parses the MSA annotations to return a Dict from the sequence identifier in the MSA to PDB and chain codes.  getseq2pdb(msa)Once you know the association between PDB chains and sequences, you can use that information together with the msacolumn2pdbresidue function to get the PDB residue number that correspond to each MSA column for given a determined sequence and PDB chain. That function downloads information from SIFTS to generate the mapping.  col2res = msacolumn2pdbresidue(msa, \"MAA_ECOLI/7-58\", \"1OCX\", \"C\")The returned dictionary can be used to get the PDB residue associated to each column (using the msaresidues function)...  using MIToS.PDB\npdbfile = downloadpdb(\"1OCX\")\npdb = read(pdbfile, PDBML)\nresdict = @residuesdict pdb model \"1\" chain \"C\" group \"ATOM\" residue All\n\nmsaresidues(msa, resdict, col2res)...or to delete the columns without PDB residues (using the hasresidues function):  using MIToS.MSA\nfiltercolumns!(msa, hasresidues(msa, col2res))"
},

{
    "location": "Pfam.html#PDB-contacts-and-AUC-1",
    "page": "Pfam",
    "title": "PDB contacts and AUC",
    "category": "section",
    "text": "The Dict between MSA columns and PDB residue number also can be used to generate a protein contact map associated to the MSA.  cmap = msacontacts(msa, resdict, col2res)That protein contact map can be used to calculate the Area Under the ROC Curve for a given score with the AUC function.  using MIToS.Information\nZMIp, MIp = buslje09(msa)\n\nAUC(ZMIp, cmap)"
},

{
    "location": "Scripts.html#",
    "page": "Scripts",
    "title": "Scripts",
    "category": "page",
    "text": ""
},

{
    "location": "Scripts.html#Scripts-1",
    "page": "Scripts",
    "title": "Scripts",
    "category": "section",
    "text": "MIToS implements several useful scripts to command line execution (without requiring Julia coding). All this scripts are located in the scripts folder of the MIToS directory. You can copy them to your working directory, use the path to their folder or put them in the path (look into the Installation section of this manual).  Pages = [\"Scripts.md\"]\nDepth = 4"
},

{
    "location": "Scripts.html#Buslje09.jl-1",
    "page": "Scripts",
    "title": "Buslje09.jl",
    "category": "section",
    "text": "julia_path = ENV[\"_\"]\nscript_path = joinpath(Pkg.dir(\"MIToS\"), \"scripts\", \"Buslje09.jl\")\nrun(`$julia_path $script_path -h`)"
},

{
    "location": "Scripts.html#BLMI.jl-1",
    "page": "Scripts",
    "title": "BLMI.jl",
    "category": "section",
    "text": "julia_path = ENV[\"_\"]\nscript_path = joinpath(Pkg.dir(\"MIToS\"), \"scripts\", \"BLMI.jl\")\nrun(`$julia_path $script_path -h`)"
},

{
    "location": "Scripts.html#DownloadPDB.jl-1",
    "page": "Scripts",
    "title": "DownloadPDB.jl",
    "category": "section",
    "text": "julia_path = ENV[\"_\"]\nscript_path = joinpath(Pkg.dir(\"MIToS\"), \"scripts\", \"DownloadPDB.jl\")\nrun(`$julia_path $script_path -h`)"
},

{
    "location": "Scripts.html#Distances.jl-1",
    "page": "Scripts",
    "title": "Distances.jl",
    "category": "section",
    "text": "julia_path = ENV[\"_\"]\nscript_path = joinpath(Pkg.dir(\"MIToS\"), \"scripts\", \"Distances.jl\")\nrun(`$julia_path $script_path -h`)"
},

{
    "location": "Scripts.html#MSADescription.jl-1",
    "page": "Scripts",
    "title": "MSADescription.jl",
    "category": "section",
    "text": "julia_path = ENV[\"_\"]\nscript_path = joinpath(Pkg.dir(\"MIToS\"), \"scripts\", \"MSADescription.jl\")\nrun(`$julia_path $script_path -h`)"
},

{
    "location": "Scripts.html#PercentIdentity.jl-1",
    "page": "Scripts",
    "title": "PercentIdentity.jl",
    "category": "section",
    "text": "julia_path = ENV[\"_\"]\nscript_path = joinpath(Pkg.dir(\"MIToS\"), \"scripts\", \"PercentIdentity.jl\")\nrun(`$julia_path $script_path -h`)"
},

{
    "location": "Scripts.html#AlignedColumns.jl-1",
    "page": "Scripts",
    "title": "AlignedColumns.jl",
    "category": "section",
    "text": "julia_path = ENV[\"_\"]\nscript_path = joinpath(Pkg.dir(\"MIToS\"), \"scripts\", \"AlignedColumns.jl\")\nrun(`$julia_path $script_path -h`)"
},

{
    "location": "Scripts.html#SplitStockholm.jl-1",
    "page": "Scripts",
    "title": "SplitStockholm.jl",
    "category": "section",
    "text": "julia_path = ENV[\"_\"]\nscript_path = joinpath(Pkg.dir(\"MIToS\"), \"scripts\", \"SplitStockholm.jl\")\nrun(`$julia_path $script_path -h`)"
},

{
    "location": "MSA_API.html#",
    "page": "MSA",
    "title": "MSA",
    "category": "page",
    "text": ""
},

{
    "location": "MSA_API.html#MIToS.MSA",
    "page": "MSA",
    "title": "MIToS.MSA",
    "category": "Module",
    "text": "The MSA module of MIToS has utilities for working with Multiple Sequence Alignments of protein Sequences (MSA).\n\nFeatures\n\nRead and write MSAs in Stockholm, FASTA or Raw format\nHandle MSA annotations\nEdit the MSA, e.g. delete columns or sequences, change sequence order, shuffling...\nKeep track of positions and annotations after modifications on the MSA\nDescribe a MSA, e.g. mean percent identity, sequence coverage, gap percentage...\n\nusing MIToS.MSA\n\n\n\n"
},

{
    "location": "MSA_API.html#MSA-1",
    "page": "MSA",
    "title": "MSA",
    "category": "section",
    "text": "MIToS.MSA"
},

{
    "location": "MSA_API.html#Contents-1",
    "page": "MSA",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"MSA_API.md\"]\nDepth = 2"
},

{
    "location": "MSA_API.html#MIToS.MSA.AbstractAlignedObject",
    "page": "MSA",
    "title": "MIToS.MSA.AbstractAlignedObject",
    "category": "Type",
    "text": "MIToS MSA and aligned sequences (aligned objects) are subtypes of AbstractMatrix{Residue}, because MSAs and sequences are stored as Matrix of Residues.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.AbstractAlignedSequence",
    "page": "MSA",
    "title": "MIToS.MSA.AbstractAlignedSequence",
    "category": "Type",
    "text": "A MIToS aligned sequence is an AbstractMatrix{Residue} with only 1 row/sequence.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.AbstractMultipleSequenceAlignment",
    "page": "MSA",
    "title": "MIToS.MSA.AbstractMultipleSequenceAlignment",
    "category": "Type",
    "text": "MSAs are stored as Matrix{Residue}. It's possible to use a NamedArray{Residue,2} as the most simple MSA with sequence identifiers and column names.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.AlignedSequence",
    "page": "MSA",
    "title": "MIToS.MSA.AlignedSequence",
    "category": "Type",
    "text": "An AlignedSequence wraps a NamedArray{Residue,2} with only 1 row/sequence. The NamedArray stores the sequence name and original column numbers as Strings.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.AnnotatedAlignedSequence",
    "page": "MSA",
    "title": "MIToS.MSA.AnnotatedAlignedSequence",
    "category": "Type",
    "text": "This type represent an aligned sequence, similar to AlignedSequence, but It also stores its Annotations.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.AnnotatedMultipleSequenceAlignment",
    "page": "MSA",
    "title": "MIToS.MSA.AnnotatedMultipleSequenceAlignment",
    "category": "Type",
    "text": "This type represent an MSA, similar to MultipleSequenceAlignment, but It also stores Annotations. This annotations are used to store residue coordinates (i.e. mapping to UniProt residue numbers).\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.Annotations",
    "page": "MSA",
    "title": "MIToS.MSA.Annotations",
    "category": "Type",
    "text": "The Annotations type is basically a container for Dicts with the annotations of a multiple sequence alignment. Annotations was designed for storage of annotations of the Stockholm format.\n\nMIToS also uses MSA annotations to keep track of:\n\nModifications of the MSA (MIToS_...) as deletion of sequences or columns.\nPositions numbers in the original MSA file (column mapping: ColMap)\nPosition of the residues in the sequence (sequence mapping: SeqMap)\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.Clusters",
    "page": "MSA",
    "title": "MIToS.MSA.Clusters",
    "category": "Type",
    "text": "Data structure to represent sequence clusters. The sequence data itself is not included.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.GappedAlphabet",
    "page": "MSA",
    "title": "MIToS.MSA.GappedAlphabet",
    "category": "Type",
    "text": "This type defines the usual alphabet of the 20 natural residues and a gap character.\n\njulia> GappedAlphabet()\nMIToS.MSA.GappedAlphabet of length 21. Residues : res\"ARNDCQEGHILKMFPSTWYV-\"\n\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.MultipleSequenceAlignment",
    "page": "MSA",
    "title": "MIToS.MSA.MultipleSequenceAlignment",
    "category": "Type",
    "text": "This MSA type include a NamedArray wrapping a Matrix of Residues. The use of NamedArray allows to store sequence names and original column numbers as Strings, and fast indexing using them.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.NoClustering",
    "page": "MSA",
    "title": "MIToS.MSA.NoClustering",
    "category": "Type",
    "text": "Use NoClustering() to avoid the use of clustering where a Clusters type is needed.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.ReducedAlphabet",
    "page": "MSA",
    "title": "MIToS.MSA.ReducedAlphabet",
    "category": "Type",
    "text": "ReducedAlphabet allows the construction of reduced residue alphabets, where residues inside parenthesis belong to the same group.\n\njulia> ab = ReducedAlphabet(\"(AILMV)(RHK)(NQST)(DE)(FWY)CGP\")\nMIToS.MSA.ReducedAlphabet of length 8 : \"(AILMV)(RHK)(NQST)(DE)(FWY)CGP\"\n\njulia> ab[Residue('K')]\n2\n\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.Residue",
    "page": "MSA",
    "title": "MIToS.MSA.Residue",
    "category": "Type",
    "text": "Most of the MIToS design is created around the Residue bitstype. It has representations for the 20 natural amino acids, a value representing insertions and deletions (GAP, '-') and one representing unknown, ambiguous and non standard residues (XAA, 'X'). Each Residue is encoded as an integer number, with the same bit representation and size than a Int. This allows fast indexing operation of probability or frequency matrices.\n\nResidue creation and conversion\n\nCreation and conversion of Residues should be treated carefully. Residue is encoded as a 32 or 64 bits type similar to Int, to get fast indexing using Int(x::Residue). Int simply calls reinterpret without checking if the residue is valid. Valid residues have integer values in the closed interval [1,22]. convert from Int  and Char always returns valid residues, however it's possible to find invalid residues (they are shown using the character '�') after the creation of uninitialized Residue arrays (i.e. using Array). You can use zeros, ones or rand to get initialized Residue arrays with valid residues. Conversions to and from Chars changes the bit representation and allows the use of the usual character representation of residues and amino acids. This conversions are used in IO operations and always return valid residues. In conversions from Char, lowercase letters, '*', '-' and '.' are translated to GAP, letters representing the 20 natural amino (ARNDCQEGHILKMFPSTWYV) acids are translated to their corresponding Residue and any other character is translated to XAA. Since lowercase letters and dots are translated to gaps, Pfam MSA insert columns are converted to columns full of gaps.\n\njulia> alanine = Residue('A')\nA\n\njulia> Char(alanine)\n'A'\n\njulia> for residue in res\"ARNDCQEGHILKMFPSTWYV-X\"\n           println(residue, \" \", Int(residue))\n       end\nA 1\nR 2\nN 3\nD 4\nC 5\nQ 6\nE 7\nG 8\nH 9\nI 10\nL 11\nK 12\nM 13\nF 14\nP 15\nS 16\nT 17\nW 18\nY 19\nV 20\n- 21\nX 22\n\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.ResidueAlphabet",
    "page": "MSA",
    "title": "MIToS.MSA.ResidueAlphabet",
    "category": "Type",
    "text": "Abstract type to define residue alphabet types.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.UngappedAlphabet",
    "page": "MSA",
    "title": "MIToS.MSA.UngappedAlphabet",
    "category": "Type",
    "text": "This type defines the usual alphabet of the 20 natural residues, without the gap character.\n\njulia> UngappedAlphabet()\nMIToS.MSA.UngappedAlphabet of length 20. Residues : res\"ARNDCQEGHILKMFPSTWYV\"\n\n\n\n\n"
},

{
    "location": "MSA_API.html#Types-1",
    "page": "MSA",
    "title": "Types",
    "category": "section",
    "text": "Modules = [MIToS.MSA]\nPrivate = false\nOrder   = [:type]"
},

{
    "location": "MSA_API.html#MIToS.MSA.GAP",
    "page": "MSA",
    "title": "MIToS.MSA.GAP",
    "category": "Constant",
    "text": "GAP is the Residue representation on MIToS for gaps ('-', insertions and deletions). Lowercase residue characters, dots and '*' are encoded as GAP in conversion from Strings and Chars. This Residue constant is encoded as Residue(21).\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.XAA",
    "page": "MSA",
    "title": "MIToS.MSA.XAA",
    "category": "Constant",
    "text": "XAA is the Residue representation for unknown, ambiguous and non standard residues. This Residue constant is encoded as Residue(22).\n\n\n\n"
},

{
    "location": "MSA_API.html#Constants-1",
    "page": "MSA",
    "title": "Constants",
    "category": "section",
    "text": "Modules = [MIToS.MSA]\nPrivate = false\nOrder   = [:constant]"
},

{
    "location": "MSA_API.html#MIToS.MSA.@res_str-Tuple{Any}",
    "page": "MSA",
    "title": "MIToS.MSA.@res_str",
    "category": "Macro",
    "text": "The MIToS macro @res_str takes a string and returns a Vector of Residues (sequence).\n\njulia> res\"MIToS\"\n5-element Array{MIToS.MSA.Residue,1}:\n M\n I\n T\n -\n S\n\n\n\n\n"
},

{
    "location": "MSA_API.html#Macros-1",
    "page": "MSA",
    "title": "Macros",
    "category": "section",
    "text": "Modules = [MIToS.MSA]\nPrivate = false\nOrder   = [:macro]"
},

{
    "location": "MSA_API.html#Base.Random.rand-Tuple{Any,Type{MIToS.MSA.Residue}}",
    "page": "MSA",
    "title": "Base.Random.rand",
    "category": "Method",
    "text": "It chooses from the 20 natural residues (it doesn't generate gaps).\n\njulia> srand(1); # Reseed the random number generator.\n\njulia> rand(Residue)\nY\n\njulia> rand(Residue, 4, 4)\n4×4 Array{MIToS.MSA.Residue,2}:\n E  L  K  P\n V  N  A  M\n I  D  A  F\n F  Y  C  P\n\n\n\n\n"
},

{
    "location": "MSA_API.html#Base.Random.shuffle!",
    "page": "MSA",
    "title": "Base.Random.shuffle!",
    "category": "Function",
    "text": "It's like Base.shuffle. When a Matrix{Residue} is used, you can indicate if the gaps should remain their positions using the last boolean argument. The previous argument should be the dimension to shuffle, 1 for shuffling residues in a sequence (row) or 2 for shuffling residues in a column.\n\njulia> msa = hcat(res\"RRE\",res\"DDK\", res\"G--\")\n3×3 Array{MIToS.MSA.Residue,2}:\n R  D  G\n R  D  -\n E  K  -\n\njulia> srand(42);\n\njulia> shuffle(msa, 1, true)\n3×3 Array{MIToS.MSA.Residue,2}:\n G  D  R\n D  R  -\n E  K  -\n\njulia> srand(42);\n\njulia> shuffle(msa, 1, false)\n3×3 Array{MIToS.MSA.Residue,2}:\n G  D  R\n D  -  R\n -  E  K\n\n\n\n\n"
},

{
    "location": "MSA_API.html#Base.Random.shuffle-Tuple{AbstractRNG,Array{MIToS.MSA.Residue,2},Vararg{Any,N}}",
    "page": "MSA",
    "title": "Base.Random.shuffle",
    "category": "Method",
    "text": "It's like shuffle but in-place. When a Matrix{Residue} or a AbstractAlignedObject (sequence or MSA) is used, you can indicate if the gaps should remain their positions using the last boolean argument.\n\n\n\n"
},

{
    "location": "MSA_API.html#Base.isvalid-Tuple{Type{MIToS.MSA.Residue},MIToS.MSA.Residue}",
    "page": "MSA",
    "title": "Base.isvalid",
    "category": "Method",
    "text": "isvalid(res::Residue)\n\nIt returns true if the encoded integer is in the closed interval [1,22].\n\n\n\n"
},

{
    "location": "MSA_API.html#Base.names-Tuple{MIToS.MSA.ReducedAlphabet}",
    "page": "MSA",
    "title": "Base.names",
    "category": "Method",
    "text": "It returns the name of each group. The name is a string with the one letter code of each residue that belong to the group.\n\njulia> ab = ReducedAlphabet(\"(AILMV)(RHK)(NQST)(DE)(FWY)CGP\")\nMIToS.MSA.ReducedAlphabet of length 8 : \"(AILMV)(RHK)(NQST)(DE)(FWY)CGP\"\n\njulia> names(ab)\n8-element Array{String,1}:\n \"AILMV\"\n \"RHK\"\n \"NQST\"\n \"DE\"\n \"FWY\"\n \"C\"    \n \"G\"\n \"P\"\n\n\n\n\n"
},

{
    "location": "MSA_API.html#Base.parse",
    "page": "MSA",
    "title": "Base.parse",
    "category": "Function",
    "text": "parse(io, format[, output; generatemapping, useidcoordinates, deletefullgaps])\n\nThe keyword argument generatemapping (false by default) indicates if the mapping of the sequences (\"SeqMap\") and columns (\"ColMap\") and the number of columns in the original MSA (\"NCol\") should be generated and saved in the annotations. If useidcoordinates is true (default: false) the sequence IDs of the form \"ID/start-end\" are parsed and used for determining the start and end positions when the mappings are generated. deletefullgaps (true by default) indicates if columns 100% gaps (generally inserts from a HMM) must be removed from the MSA.\n\n\n\n"
},

{
    "location": "MSA_API.html#Clustering.assignments-Tuple{MIToS.MSA.Clusters}",
    "page": "MSA",
    "title": "Clustering.assignments",
    "category": "Method",
    "text": "Get a vector of assignments, where the i value is the index/number of the cluster to which the i-th sequence is assigned.\n\n\n\n"
},

{
    "location": "MSA_API.html#Clustering.nclusters-Tuple{MIToS.MSA.Clusters}",
    "page": "MSA",
    "title": "Clustering.nclusters",
    "category": "Method",
    "text": "Get the number of clusters in a Clusters object.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.adjustreference",
    "page": "MSA",
    "title": "MIToS.MSA.adjustreference",
    "category": "Function",
    "text": "Creates a new matrix of residues. This function deletes positions/columns of the MSA with gaps in the reference (first) sequence.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.adjustreference!",
    "page": "MSA",
    "title": "MIToS.MSA.adjustreference!",
    "category": "Function",
    "text": "It removes positions/columns of the MSA with gaps in the reference (first) sequence.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.annotate_modification!-Tuple{MIToS.MSA.Annotations,String}",
    "page": "MSA",
    "title": "MIToS.MSA.annotate_modification!",
    "category": "Method",
    "text": "Annotates on file annotations the modifications realized by MIToS on the MSA. It always returns true, so It can be used in a boolean context.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.annotations-Tuple{MIToS.MSA.AnnotatedMultipleSequenceAlignment}",
    "page": "MSA",
    "title": "MIToS.MSA.annotations",
    "category": "Method",
    "text": "annotations returns the Annotations of an MSA or aligned sequence.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.columngapfraction-Tuple{AbstractArray{MIToS.MSA.Residue,2}}",
    "page": "MSA",
    "title": "MIToS.MSA.columngapfraction",
    "category": "Method",
    "text": "Fraction of gaps per column/position on the MSA\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.columnnames-Tuple{NamedArrays.NamedArray{MIToS.MSA.Residue,2,AT,DT}}",
    "page": "MSA",
    "title": "MIToS.MSA.columnnames",
    "category": "Method",
    "text": "columnnames(msa)\n\nIt returns a Vector{String} with the sequence names/identifiers. If the msa is a Matrix{Residue} this function returns the actual column numbers as strings. Otherwise it returns the column number of the original MSA through the wrapped NamedArray column names.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.columnpairsmatrix-Tuple{AbstractArray{MIToS.MSA.Residue,2},Type{T},Type{Val{diagonal}},T}",
    "page": "MSA",
    "title": "MIToS.MSA.columnpairsmatrix",
    "category": "Method",
    "text": "Initialize an empty PairwiseListMatrix for a pairwise measure in sequence pairs. It uses the sequence names if they are available, otherwise it uses the actual sequence numbers. You can use the positional argument to indicate the number Type (default: Float64), if the PairwiseListMatrix should store the diagonal values on the list (default: false) and a default value for the diagonal (default: NaN).\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.coverage-Tuple{AbstractArray{MIToS.MSA.Residue,2}}",
    "page": "MSA",
    "title": "MIToS.MSA.coverage",
    "category": "Method",
    "text": "Coverage of the sequences with respect of the number of positions on the MSA\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.delete_annotated_modifications!-Tuple{MIToS.MSA.Annotations}",
    "page": "MSA",
    "title": "MIToS.MSA.delete_annotated_modifications!",
    "category": "Method",
    "text": "Deletes all the MIToS annotated modifications\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.deletefullgapcolumns!",
    "page": "MSA",
    "title": "MIToS.MSA.deletefullgapcolumns!",
    "category": "Function",
    "text": "Deletes columns with 100% gaps, this columns are generated by inserts.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.filtercolumns!",
    "page": "MSA",
    "title": "MIToS.MSA.filtercolumns!",
    "category": "Function",
    "text": "filtercolumns!(msa, mask[, annotate::Bool=true])\n\nIt allows to filter MSA or aligned sequence columns/positions using a AbstractVector{Bool} mask. Annotations are updated if annotate is true (default).\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.filtercolumns!-Tuple{MIToS.MSA.Annotations,Any}",
    "page": "MSA",
    "title": "MIToS.MSA.filtercolumns!",
    "category": "Method",
    "text": "filtercolumns!(data::Annotations, mask)\n\nIt is useful for deleting column annotations (creating a subset in place).\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.filtercolumns-Tuple{AbstractArray{MIToS.MSA.Residue,2},Any}",
    "page": "MSA",
    "title": "MIToS.MSA.filtercolumns",
    "category": "Method",
    "text": "It's similar to filtercolumns! but for an AbstractMatrix{Residue}\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.filtersequences!",
    "page": "MSA",
    "title": "MIToS.MSA.filtersequences!",
    "category": "Function",
    "text": "filtersequences!(msa, mask[, annotate::Bool=true])\n\nIt allows to filter msa sequences using a AbstractVector{Bool} mask (It removes sequnces with false values). AnnotatedMultipleSequenceAlignment annotations are updated if annotate is true (default).\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.filtersequences!-Tuple{MIToS.MSA.Annotations,Array{String,1},AbstractArray{Bool,1}}",
    "page": "MSA",
    "title": "MIToS.MSA.filtersequences!",
    "category": "Method",
    "text": "filtersequences!(data::Annotations, ids::Vector{String}, mask::AbstractArray{Bool,1})\n\nIt is useful for deleting sequence annotations. ids should be a list of the sequence names and mask should be a logical vector.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.filtersequences-Tuple{AbstractArray{MIToS.MSA.Residue,2},Any}",
    "page": "MSA",
    "title": "MIToS.MSA.filtersequences",
    "category": "Method",
    "text": "It's similar to filtersequences! but for an AbstractMatrix{Residue}\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.gapfraction-Tuple{AbstractArray{MIToS.MSA.Residue,N}}",
    "page": "MSA",
    "title": "MIToS.MSA.gapfraction",
    "category": "Method",
    "text": "It calculates the fraction of gaps on the Array (alignment, sequence, column, etc.). This function can take an extra dimension argument for calculation of the gap fraction over the given dimension.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.gapstrip!",
    "page": "MSA",
    "title": "MIToS.MSA.gapstrip!",
    "category": "Function",
    "text": "This functions deletes/filters sequences and columns/positions on the MSA on the following order:\n\nRemoves all the columns/position on the MSA with gaps on the reference (first) sequence.\nRemoves all the sequences with a coverage with respect to the number of\n\ncolumns/positions on the MSA less than a coveragelimit (default to 0.75: sequences with 25% of gaps).\n\nRemoves all the columns/position on the MSA with more than a gaplimit\n\n(default to 0.5: 50% of gaps).\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.gapstrip-Tuple{AbstractArray{MIToS.MSA.Residue,2}}",
    "page": "MSA",
    "title": "MIToS.MSA.gapstrip",
    "category": "Method",
    "text": "Creates a new matrix of Residues (MSA) with deleted sequences and columns/positions. The MSA is edited in the following way:\n\nRemoves all the columns/position on the MSA with gaps on the reference (first) sequence\nRemoves all the sequences with a coverage with respect to the number of\n\ncolumns/positions on the MSA less than a coveragelimit  (default to 0.75: sequences with 25% of gaps)\n\nRemoves all the columns/position on the MSA with more than a gaplimit\n\n(default to 0.5: 50% of gaps)\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getannotcolumn",
    "page": "MSA",
    "title": "MIToS.MSA.getannotcolumn",
    "category": "Function",
    "text": "getannotcolumn(ann[, feature[,default]])\n\nIt returns per column annotation for feature\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getannotfile",
    "page": "MSA",
    "title": "MIToS.MSA.getannotfile",
    "category": "Function",
    "text": "getannotfile(ann[, feature[,default]])\n\nIt returns per file annotation for feature\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getannotresidue",
    "page": "MSA",
    "title": "MIToS.MSA.getannotresidue",
    "category": "Function",
    "text": "getannotresidue(ann[, seqname, feature[,default]])\n\nIt returns per residue annotation for (seqname, feature)\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getannotsequence",
    "page": "MSA",
    "title": "MIToS.MSA.getannotsequence",
    "category": "Function",
    "text": "getannotsequence(ann[, seqname, feature[,default]])\n\nIt returns per sequence annotation for (seqname, feature)\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getcolumnmapping-Tuple{MIToS.MSA.AnnotatedMultipleSequenceAlignment}",
    "page": "MSA",
    "title": "MIToS.MSA.getcolumnmapping",
    "category": "Method",
    "text": "It returns a Vector{Int} with the original column number of each column on the actual MSA. The mapping is annotated in the \"ColMap\" file annotation of an AnnotatedMultipleSequenceAlignment or in the column names of an NamedArray or MultipleSequenceAlignment.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getnamedict-Tuple{MIToS.MSA.ReducedAlphabet}",
    "page": "MSA",
    "title": "MIToS.MSA.getnamedict",
    "category": "Method",
    "text": "It takes a ResidueAlphabet and returns a dictionary from group name to group position.\n\njulia> ab = ReducedAlphabet(\"(AILMV)(RHK)(NQST)(DE)(FWY)CGP\")\nMIToS.MSA.ReducedAlphabet of length 8 : \"(AILMV)(RHK)(NQST)(DE)(FWY)CGP\"\n\njulia> getnamedict(ab)\nDataStructures.OrderedDict{String,Int64} with 8 entries:\n  \"AILMV\" => 1\n  \"RHK\"   => 2\n  \"NQST\"  => 3\n  \"DE\"    => 4\n  \"FWY\"   => 5\n  \"C\"     => 6\n  \"G\"     => 7\n  \"P\"     => 8\n\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getresidues-Tuple{Array{MIToS.MSA.Residue,2}}",
    "page": "MSA",
    "title": "MIToS.MSA.getresidues",
    "category": "Method",
    "text": "getresidues allows you to access the residues stored inside an MSA or aligned sequence as a Matrix{Residue} without annotations nor column/row names.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getresiduesequences-Tuple{Array{MIToS.MSA.Residue,2}}",
    "page": "MSA",
    "title": "MIToS.MSA.getresiduesequences",
    "category": "Method",
    "text": "getresiduesequences returns a Vector{Vector{Residue}} with all the MSA sequences without annotations nor column/sequence names.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getsequence",
    "page": "MSA",
    "title": "MIToS.MSA.getsequence",
    "category": "Function",
    "text": "getsequence takes an MSA and a sequence number or identifier and returns an aligned sequence object. If the MSA is an AnnotatedMultipleSequenceAlignment, it returns an AnnotatedAlignedSequence with the sequence annotations. From a MultipleSequenceAlignment, It returns an AlignedSequence object. If an Annotations object and a sequence identifier are used, this function returns the annotations related to the sequence.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getsequencemapping-Tuple{MIToS.MSA.AnnotatedMultipleSequenceAlignment,String}",
    "page": "MSA",
    "title": "MIToS.MSA.getsequencemapping",
    "category": "Method",
    "text": "It returns the sequence coordinates as a Vector{Int} for an MSA sequence. That vector has one element for each MSA column. If the number if 0 in the mapping, there is a gap in that column for that sequence.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.getweight-Tuple{MIToS.MSA.NoClustering,Int64}",
    "page": "MSA",
    "title": "MIToS.MSA.getweight",
    "category": "Method",
    "text": "getweight(c[, i::Int])\n\nThis function returns the weight of the sequence number i. getweight should be defined for any type used for count!/count in order to use his weigths. If i isn't used, this function returns a vector with the weight of each sequence.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.hobohmI-Tuple{AbstractArray{MIToS.MSA.Residue,2},Any}",
    "page": "MSA",
    "title": "MIToS.MSA.hobohmI",
    "category": "Method",
    "text": "Sequence clustering using the Hobohm I method from Hobohm et. al. 1992.\n\nHobohm, Uwe, et al. \"Selection of representative protein data sets.\" Protein Science 1.3 (1992): 409-417.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.meanpercentidentity",
    "page": "MSA",
    "title": "MIToS.MSA.meanpercentidentity",
    "category": "Function",
    "text": "Returns the mean of the percent identity between the sequences of a MSA. If the MSA has 300 sequences or less, the mean is exact. If the MSA has more sequences and the exact keyword is false (defualt), 44850 random pairs of sequences are used for the estimation. The number of samples can be changed using the second argument. Use exact=true to perform all the pairwise comparison (the calculation could be slow).\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.namedmatrix-Tuple{MIToS.MSA.AbstractAlignedObject}",
    "page": "MSA",
    "title": "MIToS.MSA.namedmatrix",
    "category": "Method",
    "text": "namedmatrix returns the NamedArray{Residue,2} stored in an MSA or aligned sequence.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.ncolumns-Tuple{AbstractArray{MIToS.MSA.Residue,2}}",
    "page": "MSA",
    "title": "MIToS.MSA.ncolumns",
    "category": "Method",
    "text": "ncolumns returns the number of MSA columns or positions.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.ncolumns-Tuple{MIToS.MSA.Annotations}",
    "page": "MSA",
    "title": "MIToS.MSA.ncolumns",
    "category": "Method",
    "text": "ncolumns(ann::Annotations) returns the number of columns/residues with annotations. This function returns -1 if there is not annotations per column/residue.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.nsequences-Tuple{AbstractArray{MIToS.MSA.Residue,2}}",
    "page": "MSA",
    "title": "MIToS.MSA.nsequences",
    "category": "Method",
    "text": "nsequences returns the number of sequences on the MSA.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.percentidentity",
    "page": "MSA",
    "title": "MIToS.MSA.percentidentity",
    "category": "Function",
    "text": "percentidentity(msa[, out::Type=Float64])\n\nCalculates the identity between all the sequences on a MSA. You can indicate the output element type with the last optional parameter (Float64 by default). For a MSA with a lot of sequences, you can use Float32 or Flot16 in order to avoid the OutOfMemoryError().\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.percentidentity-Tuple{Any,Any,Any}",
    "page": "MSA",
    "title": "MIToS.MSA.percentidentity",
    "category": "Method",
    "text": "percentidentity(seq1, seq2, threshold)\n\nComputes quickly if two aligned sequences have a identity value greater than a given threshold value. Returns a boolean value. Positions with gaps in both sequences doesn't count to the length of the sequences. Positions with a XAA in at least one sequence aren't counted.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.percentidentity-Tuple{Any,Any}",
    "page": "MSA",
    "title": "MIToS.MSA.percentidentity",
    "category": "Method",
    "text": "percentidentity(seq1, seq2)\n\nCalculates the fraction of identities between two aligned sequences. The identity value is calculated as the number of identical characters in the i-th position of both sequences divided by the length of both sequences. Positions with gaps in both sequences doesn't count to the length of the sequences. Positions with a XAA in at least one sequence aren't counted.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.percentsimilarity",
    "page": "MSA",
    "title": "MIToS.MSA.percentsimilarity",
    "category": "Function",
    "text": "Calculates the similarity percent between two aligned sequences. The 100% is the length of the aligned sequences minus the number of columns with gaps in both sequences and the number of columns with at least one residue outside the alphabet. So, columns with residues outside the alphabet (other than the specially treated GAP) aren't counted to the protein length. Two residues are considered similar if they below to the same group in a ReducedAlphabet. The alphabet (third positional argument) by default is:\n\nReducedAlphabet(\"(AILMV)(NQST)(RHK)(DE)(FWY)CGP\")\n\nThe first group is composed of the non polar residues (AILMV), the second group is composed of polar residues, the third group are positive residues, the fourth group are negative residues, the fifth group is composed by the aromatic residues (FWY). C, G and P are considered unique residues.\n\nOther residue groups/alphabets:\n\nSMS (Sequence Manipulation Suite) Ident and Sim:\n\nReducedAlphabet(\"(GAVLI)(FYW)(ST)(KRH)(DENQ)P(CM)\")\n\nStothard P (2000) The Sequence Manipulation Suite: JavaScript programs for analyzing and formatting protein and DNA sequences. Biotechniques 28:1102-1104.\n\nBio3D 2.2 seqidentity:\n\nReducedAlphabet(\"(GA)(MVLI)(FYW)(ST)(KRH)(DE)(NQ)PC\")\n\nGrant, B.J. et al. (2006) Bioinformatics 22, 2695–2696.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.percentsimilarity-Tuple{AbstractArray{MIToS.MSA.Residue,2},Vararg{Any,N}}",
    "page": "MSA",
    "title": "MIToS.MSA.percentsimilarity",
    "category": "Method",
    "text": "Calculates the similarity percent between all the sequences on a MSA. You can indicate the output element type with the out keyword argument (Float64 by default). For an MSA with a lot of sequences, you can use out=Float32 or out=Flot16 in order to avoid the OutOfMemoryError().\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.printmodifications-Tuple{MIToS.MSA.Annotations}",
    "page": "MSA",
    "title": "MIToS.MSA.printmodifications",
    "category": "Method",
    "text": "Prints MIToS annotated modifications\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.residue2three-Tuple{MIToS.MSA.Residue}",
    "page": "MSA",
    "title": "MIToS.MSA.residue2three",
    "category": "Method",
    "text": "This function returns the three letter name of the Residue.\n\njulia> residue2three(Residue('G'))\n\"GLY\"\n\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.residuefraction-Tuple{AbstractArray{MIToS.MSA.Residue,N}}",
    "page": "MSA",
    "title": "MIToS.MSA.residuefraction",
    "category": "Method",
    "text": "It calculates the fraction of residues (no gaps) on the Array (alignment, sequence, column, etc.). This function can take an extra dimension argument for calculation of the residue fraction over the given dimension\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.sequencenames-Tuple{NamedArrays.NamedArray{MIToS.MSA.Residue,2,AT,DT}}",
    "page": "MSA",
    "title": "MIToS.MSA.sequencenames",
    "category": "Method",
    "text": "sequencenames(msa)\n\nIt returns a Vector{String} with the sequence names/identifiers.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.sequencepairsmatrix-Tuple{AbstractArray{MIToS.MSA.Residue,2},Type{T},Type{Val{diagonal}},T}",
    "page": "MSA",
    "title": "MIToS.MSA.sequencepairsmatrix",
    "category": "Method",
    "text": "Initialize an empty PairwiseListMatrix for a pairwise measure in column pairs. It uses the column mapping (column number in the input MSA file) if it’s available, otherwise it uses the actual column numbers. You can use the positional argument to indicate the number Type (default: Float64), if the PairwiseListMatrix should store the diagonal values on the list (default: false) and a default value for the diagonal (default: NaN).\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.setannotcolumn!",
    "page": "MSA",
    "title": "MIToS.MSA.setannotcolumn!",
    "category": "Function",
    "text": "setannotcolumn!(ann, feature, annotation)\n\nIt stores per column annotation (1 char per column) for feature\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.setannotfile!",
    "page": "MSA",
    "title": "MIToS.MSA.setannotfile!",
    "category": "Function",
    "text": "setannotfile!(ann, feature, annotation)\n\nIt stores per file annotation for feature\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.setannotresidue!",
    "page": "MSA",
    "title": "MIToS.MSA.setannotresidue!",
    "category": "Function",
    "text": "setannotresidue!(ann, seqname, feature, annotation)\n\nIt stores per residue annotation (1 char per residue) for (seqname, feature)\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.setannotsequence!",
    "page": "MSA",
    "title": "MIToS.MSA.setannotsequence!",
    "category": "Function",
    "text": "setannotsequence!(ann, seqname, feature, annotation)\n\nIt stores per sequence annotation for (seqname, feature)\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.setreference!",
    "page": "MSA",
    "title": "MIToS.MSA.setreference!",
    "category": "Function",
    "text": "It puts the sequence i (name or position) as reference (first sequence) of the MSA. This function swaps the sequences 1 and i.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.stringsequence-Tuple{AbstractArray{MIToS.MSA.Residue,2},Any}",
    "page": "MSA",
    "title": "MIToS.MSA.stringsequence",
    "category": "Method",
    "text": "stringsequence(seq)\nstringsequence(msa, i::Int)\nstringsequence(msa, id::String)\n\nIt returns the selected sequence as a String.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.swapsequences!-Tuple{Array{MIToS.MSA.Residue,2},Int64,Int64}",
    "page": "MSA",
    "title": "MIToS.MSA.swapsequences!",
    "category": "Method",
    "text": "It swaps the sequences on the positions i and j of an MSA. Also it's possible to swap sequences using their sequence names/identifiers when the MSA object as names.\n\n\n\n"
},

{
    "location": "MSA_API.html#MIToS.MSA.three2residue-Tuple{String}",
    "page": "MSA",
    "title": "MIToS.MSA.three2residue",
    "category": "Method",
    "text": "It takes a three letter residue name and returns the corresponding Residue. If the name isn't in the MIToS dictionary, a XAA is returned.\n\njulia> three2residue(\"ALA\")\nA\n\n\n\n\n"
},

{
    "location": "MSA_API.html#StatsBase.counts-Tuple{MIToS.MSA.Clusters}",
    "page": "MSA",
    "title": "StatsBase.counts",
    "category": "Method",
    "text": "Get sample counts of clusters as a Vector. Each k value is the number of samples assigned to the k-th cluster.\n\n\n\n"
},

{
    "location": "MSA_API.html#Methods-and-functions-1",
    "page": "MSA",
    "title": "Methods and functions",
    "category": "section",
    "text": "Modules = [MIToS.MSA]\nPrivate = false\nOrder   = [:function]"
},

{
    "location": "Information_API.html#",
    "page": "Information",
    "title": "Information",
    "category": "page",
    "text": ""
},

{
    "location": "Information_API.html#MIToS.Information",
    "page": "Information",
    "title": "MIToS.Information",
    "category": "Module",
    "text": "The Information module of MIToS defines types and functions useful to calculate information measures (e.g. Mutual Information (MI) and Entropy) over a Multiple Sequence Alignment (MSA). This module was designed to count Residues (defined in the MSA module) in special contingency tables (as fast as possible) and to derive probabilities from this counts. Also, includes methods for applying corrections to that tables, e.g. pseudocounts and pseudo frequencies. Finally, Information allows to use this probabilities and counts to estimate information measures and other frequency based values.\n\nFeatures\n\nEstimate multi dimensional frequencies and probabilities tables from sequences, MSAs, etc...\nCorrection for small number of observations\nCorrection for data redundancy on a MSA\nEstimate information measures\nCalculate corrected mutual information between residues\n\nusing MIToS.Information\n\n\n\n"
},

{
    "location": "Information_API.html#Information-1",
    "page": "Information",
    "title": "Information",
    "category": "section",
    "text": "MIToS.Information"
},

{
    "location": "Information_API.html#Contents-1",
    "page": "Information",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"Information_API.md\"]\nDepth = 2"
},

{
    "location": "Information_API.html#MIToS.Information.AdditiveSmoothing",
    "page": "Information",
    "title": "MIToS.Information.AdditiveSmoothing",
    "category": "Type",
    "text": "Additive Smoothing or fixed pseudocount λ for ResidueCount (in order to estimate probabilities when the number of samples is low).\n\nCommon values of λ are:\n\n0 : No cell frequency prior, gives you the maximum likelihood estimator.\n0.05 is the optimum value for λ found in Buslje et. al. 2009, similar results was obtained for λ in the range [0.025, 0.075].\n1 / p : Perks prior (Perks, 1947) where p the number of parameters (i.e. residues, pairs of residues) to estimate. If p is the number of residues (20 without counting gaps), this gives you 0.05.\nsqrt(n) / p : Minimax prior (Trybula, 1958) where n is the number of samples and p the number of parameters to estimate. If the number of samples n is 400 (minimum number of sequence clusters for achieve good performance in Buslje et. al. 2009) for estimating 400 parameters (pairs of residues without counting gaps) this gives you 0.05.\n0.5 : Jeffreys prior (Jeffreys, 1946).\n1 : Bayes-Laplace uniform prior, aka. Laplace smoothing.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.BLOSUM_Pseudofrequencies",
    "page": "Information",
    "title": "MIToS.Information.BLOSUM_Pseudofrequencies",
    "category": "Type",
    "text": "BLOSUM_Pseudofrequencies type. It takes to arguments/fields:\n\nα : Usually the number of sequences or sequence clusters in the MSA.\nβ : The weight of the pseudofrequencies, a value close to 8.512 when α is the number of sequence clusters.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.ContingencyTable",
    "page": "Information",
    "title": "MIToS.Information.ContingencyTable",
    "category": "Type",
    "text": "A ContingencyTable is a multidimensional array. It stores the contingency matrix, its marginal values and total. The type also has an internal and private temporal array and an alphabet object. It's a parametric type, taking three ordered parameters:\n\nT : The element type of the multidimensional array.\nN : It's the dimension of the array and should be an Int.\nA : This should be a type, subtype of ResidueAlphabet, i.e.: UngappedAlphabet,\n\nGappedAlphabet or ReducedAlphabet.\n\nA ContingencyTable can be created from an alphabet if all the parameters are given. Otherwise, you need to give a type, a number (Val) and an alphabet. You can also create a ContingencyTable using a matrix and a alphabet. For example:\n\nContingencyTable{Float64, 2, UngappedAlphabet}(UngappedAlphabet())\nContingencyTable(Float64, Val{2}, UngappedAlphabet())\nContingencyTable(zeros(Float64,20,20), UngappedAlphabet())\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.Counts",
    "page": "Information",
    "title": "MIToS.Information.Counts",
    "category": "Type",
    "text": "A Counts object wraps a ContingencyTable storing counts/frequencies.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.NoPseudocount",
    "page": "Information",
    "title": "MIToS.Information.NoPseudocount",
    "category": "Type",
    "text": "You can use NoPseudocount() to avoid pseudocount corrections where a Pseudocount type is needed.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.NoPseudofrequencies",
    "page": "Information",
    "title": "MIToS.Information.NoPseudofrequencies",
    "category": "Type",
    "text": "You can use NoPseudofrequencies() to avoid pseudocount corrections where a Pseudofrequencies type is needed.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.Probabilities",
    "page": "Information",
    "title": "MIToS.Information.Probabilities",
    "category": "Type",
    "text": "A Probabilities object wraps a ContingencyTable storing probabilities. It doesn't perform any check. If the total isn't one, you must use normalize or normalize!on the ContingencyTable before wrapping it to make the sum of the probabilities equal to one.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.Pseudocount",
    "page": "Information",
    "title": "MIToS.Information.Pseudocount",
    "category": "Type",
    "text": "Parametric abstract type to define pseudocount types\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.Pseudofrequencies",
    "page": "Information",
    "title": "MIToS.Information.Pseudofrequencies",
    "category": "Type",
    "text": "Parametric abstract type to define pseudofrequencies types\n\n\n\n"
},

{
    "location": "Information_API.html#Types-1",
    "page": "Information",
    "title": "Types",
    "category": "section",
    "text": "Modules = [MIToS.Information]\nPrivate = false\nOrder   = [:type]"
},

{
    "location": "Information_API.html#MIToS.Information.BLOSUM62_Pi",
    "page": "Information",
    "title": "MIToS.Information.BLOSUM62_Pi",
    "category": "Constant",
    "text": "BLOSUM62 probabilities P(aa) for each residue on the UngappedAlphabet. SUM:  0.9987\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.BLOSUM62_Pij",
    "page": "Information",
    "title": "MIToS.Information.BLOSUM62_Pij",
    "category": "Constant",
    "text": "Table with conditional probabilities of residues based on BLOSUM62. The normalization is done row based. The firts row contains the P(aa|A) and so one.\n\n\n\n"
},

{
    "location": "Information_API.html#Constants-1",
    "page": "Information",
    "title": "Constants",
    "category": "section",
    "text": "Modules = [MIToS.Information]\nPrivate = false\nOrder   = [:constant]"
},

{
    "location": "Information_API.html#Macros-1",
    "page": "Information",
    "title": "Macros",
    "category": "section",
    "text": "Modules = [MIToS.Information]\nPrivate = false\nOrder   = [:macro]"
},

{
    "location": "Information_API.html#Base.LinAlg.normalize!-Tuple{MIToS.Information.ContingencyTable{T,N,A}}",
    "page": "Information",
    "title": "Base.LinAlg.normalize!",
    "category": "Method",
    "text": "normalize! makes the sum of the frequencies to be one, in place.\n\n\n\n"
},

{
    "location": "Information_API.html#Base.LinAlg.normalize-Tuple{MIToS.Information.ContingencyTable{T,N,A}}",
    "page": "Information",
    "title": "Base.LinAlg.normalize",
    "category": "Method",
    "text": "normalize returns another table where the sum of the frequencies is one.\n\n\n\n"
},

{
    "location": "Information_API.html#Base.count-Tuple{Vararg{AbstractArray{MIToS.MSA.Residue,1},N}}",
    "page": "Information",
    "title": "Base.count",
    "category": "Method",
    "text": "It returns a ContingencyTable wrapped in a Counts type with the frequencies of residues in the sequences that takes as arguments. The dimension of the table is equal to the number of sequences. You can use the keyword arguments alphabet, weights and pseudocounts to indicate the alphabet of the table (default to UngappedAlphabet()), a clustering result (default to NoClustering()) and the pseudocounts (default to NoPseudocount()) to be used during the estimation of the frequencies.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.APC!-Tuple{Array{T,2}}",
    "page": "Information",
    "title": "MIToS.Information.APC!",
    "category": "Method",
    "text": "APC (Dunn et. al. 2008)\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.BLMI-Tuple{AbstractArray{MIToS.MSA.Residue,2}}",
    "page": "Information",
    "title": "MIToS.Information.BLMI",
    "category": "Method",
    "text": "BLMI takes a MSA or a file and a Format as first arguments. It calculates a Z score (ZBLMI) and a corrected MI/MIp as described on Busjle et. al. 2009 but using using BLOSUM62 pseudo frequencies instead of a fixed pseudocount.\n\nKeyword argument, type, default value and descriptions:\n\n  - beta        Float64   8.512   β for BLOSUM62 pseudo frequencies\n  - lambda      Float64   0.0     Low count value\n  - threshold             62      Percent identity threshold for sequence clustering (Hobohm I)\n  - maxgap      Float64   0.5     Maximum fraction of gaps in positions included in calculation\n  - apc         Bool      true    Use APC correction (MIp)\n  - samples     Int       50      Number of samples for Z-score\n  - fixedgaps   Bool      true    Fix gaps positions for the random samples\n\nThis function returns:\n\n  - Z score (ZBLMI)\n  - MI or MIp using BLOSUM62 pseudo frequencies (BLMI/BLMIp)\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.apply_pseudocount!-Tuple{MIToS.Information.ContingencyTable{T,N,A},T}",
    "page": "Information",
    "title": "MIToS.Information.apply_pseudocount!",
    "category": "Method",
    "text": "It adds the pseudocount value to the table cells.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.apply_pseudofrequencies!-Tuple{MIToS.Information.ContingencyTable{T,2,MIToS.MSA.UngappedAlphabet},MIToS.Information.BLOSUM_Pseudofrequencies}",
    "page": "Information",
    "title": "MIToS.Information.apply_pseudofrequencies!",
    "category": "Method",
    "text": "apply_pseudofrequencies!{T}(Pab::ContingencyTable{T,2,UngappedAlphabet}, pseudofrequencies::BLOSUM_Pseudofrequencies)\n\nWhen a BLOSUM_Pseudofrequencies(α,β) is used, this function applies pseudofrequencies Gab over Pab, as a weighted mean of both. It uses the conditional probability matrix BLOSUM62_Pij and the real frequencies/probabilities Pab to estimate the pseudofrequencies Gab. α is the weight of the real frequencies Pab and β the weight of the pseudofrequencies.\n\nGab = Σcd  Pcd ⋅ BLOSUM62( a | c ) ⋅ BLOSUM62( b | d ) Pab = (α ⋅ Pab + β ⋅ Gab )/(α + β)\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.buslje09-Tuple{AbstractArray{MIToS.MSA.Residue,2}}",
    "page": "Information",
    "title": "MIToS.Information.buslje09",
    "category": "Method",
    "text": "buslje09 takes a MSA or a file and a Format as first arguments. It calculates a Z score and a corrected MI/MIp as described on Busjle et. al. 2009.\n\nkeyword argument, type, default value and descriptions:\n\n  - lambda      Float64   0.05    Low count value\n  - clustering  Bool      true    Sequence clustering (Hobohm I)\n  - threshold             62      Percent identity threshold for clustering\n  - maxgap      Float64   0.5     Maximum fraction of gaps in positions included in calculation\n  - apc         Bool      true    Use APC correction (MIp)\n  - samples     Int       100     Number of samples for Z-score\n  - fixedgaps   Bool      true    Fix gaps positions for the random samples\n  - alphabet    ResidueAlphabet UngappedAlphabet()  Residue alphabet to be used\n\nThis function returns:\n\n  - Z score\n  - MI or MIp\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.count!-Tuple{MIToS.Information.ContingencyTable{T,N,A},Any,MIToS.Information.Pseudocount,Vararg{AbstractArray{MIToS.MSA.Residue,1},N}}",
    "page": "Information",
    "title": "MIToS.Information.count!",
    "category": "Method",
    "text": "It populates a ContingencyTable (first argument) using the frequencies in the sequences (last positional arguments). The dimension of the table must match the number of sequences and all the sequences must have the same length. You must indicate the used weights and pseudocounts as second and third positional arguments respectively. You can use NoPseudofrequencies() and NoClustering() to avoid the use of sequence weighting and pseudocounts, respectively.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.cumulative-Tuple{PairwiseListMatrices.PairwiseListMatrix{T,D,VT},T}",
    "page": "Information",
    "title": "MIToS.Information.cumulative",
    "category": "Method",
    "text": "cumulative allows to calculate cumulative scores (i.e. cMI) as defined in Buslje et. al. 2010\n\n\"We calculated a cumulative mutual information score (cMI) for each residue as the sum of MI values above a certain threshold for every amino acid pair where the particular residue appears. This value defines to what degree a given amino acid takes part in a mutual information network.\" Buslje, Cristina Marino, Elin Teppa, Tomas Di Doménico, José María Delfino, and Morten Nielsen. Networks of high mutual information define the structural proximity of catalytic sites: implications for catalytic residue identification. PLoS Comput Biol 6, no. 11 (2010): e1000978.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.delete_dimensions!-Tuple{MIToS.Information.ContingencyTable{T,S,A},MIToS.Information.ContingencyTable{T,N,A},Vararg{Int64,N}}",
    "page": "Information",
    "title": "MIToS.Information.delete_dimensions!",
    "category": "Method",
    "text": "delete_dimensions!(out::ContingencyTable, in::ContingencyTable, dimensions::Int...)\n\nThis function fills a ContingencyTable with the counts/probabilities on in after the deletion of dimensions. i.e. This is useful for getting Pxy from Pxyz.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.delete_dimensions-Tuple{MIToS.Information.ContingencyTable{T,N,A},Vararg{Int64,I}}",
    "page": "Information",
    "title": "MIToS.Information.delete_dimensions",
    "category": "Method",
    "text": "delete_dimensions(in::ContingencyTable, dimensions::Int...)\n\nThis function creates a ContingencyTable with the counts/probabilities on in after the deletion of dimensions. i.e. This is useful for getting Pxy from Pxyz.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.gap_intersection_percentage-Tuple{MIToS.Information.Counts{T,2,MIToS.MSA.GappedAlphabet}}",
    "page": "Information",
    "title": "MIToS.Information.gap_intersection_percentage",
    "category": "Method",
    "text": "It calculates the gap intersection as percentage from a table of Counts.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.gap_union_percentage-Tuple{MIToS.Information.Counts{T,2,MIToS.MSA.GappedAlphabet}}",
    "page": "Information",
    "title": "MIToS.Information.gap_union_percentage",
    "category": "Method",
    "text": "It calculates the gap union as percentage from a table of Counts.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.gaussdca-Tuple{Any}",
    "page": "Information",
    "title": "MIToS.Information.gaussdca",
    "category": "Method",
    "text": "Wrapper function to GaussDCA.gDCA. You need to install GaussDCA:\n\nPkg.clone(\"https://github.com/carlobaldassi/GaussDCA.jl\")\n\nLook into GaussDCA.jl README for further information. If you use this wrapper, please cite the GaussDCA publication and the package's doi.\n\nGaussDCA Publication: Baldassi, Carlo, Marco Zamparo, Christoph Feinauer, Andrea Procaccini, Riccardo Zecchina, Martin Weigt, and Andrea Pagnani. \"Fast and accurate multivariate Gaussian modeling of protein families: predicting residue contacts and protein-interaction partners.\" PloS one 9, no. 3 (2014): e92721.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.getalphabet-Tuple{MIToS.Information.ContingencyTable}",
    "page": "Information",
    "title": "MIToS.Information.getalphabet",
    "category": "Method",
    "text": "getalphabet allows to access the stored alphabet object.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.getcontingencytable-Tuple{MIToS.Information.Probabilities{T,N,A}}",
    "page": "Information",
    "title": "MIToS.Information.getcontingencytable",
    "category": "Method",
    "text": "getcontingencytable allows to access the wrapped ContingencyTable in a Probabilities or Counts object.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.getmarginals-Tuple{MIToS.Information.ContingencyTable}",
    "page": "Information",
    "title": "MIToS.Information.getmarginals",
    "category": "Method",
    "text": "getmarginals allows to access the array with the marginal values (NamedArray).\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.getmarginalsarray-Tuple{MIToS.Information.ContingencyTable}",
    "page": "Information",
    "title": "MIToS.Information.getmarginalsarray",
    "category": "Method",
    "text": "getmarginalsarray allows to access the array with the marginal values (Array without names).\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.gettable-Tuple{MIToS.Information.ContingencyTable}",
    "page": "Information",
    "title": "MIToS.Information.gettable",
    "category": "Method",
    "text": "gettable allows to access the table (NamedArray).\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.gettablearray-Tuple{MIToS.Information.ContingencyTable}",
    "page": "Information",
    "title": "MIToS.Information.gettablearray",
    "category": "Method",
    "text": "gettablearray allows to access the table (Array without names).\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.gettotal-Tuple{MIToS.Information.ContingencyTable}",
    "page": "Information",
    "title": "MIToS.Information.gettotal",
    "category": "Method",
    "text": "gettotal allows to access the stored total value.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.kullback_leibler-Tuple{MIToS.Information.Probabilities{T,1,A},Any,Real}",
    "page": "Information",
    "title": "MIToS.Information.kullback_leibler",
    "category": "Method",
    "text": "It calculates the Kullback-Leibler (KL) divergence from a table of Probabilities. The second positional argument is a Probabilities or ContingencyTable with the background distribution. It's optional, the default is the BLOSUM62_Pi table. Use last and optional positional argument to change the base of the log. The default base is e, so the result is in nats. You can use 2.0 as base to get the result in bits.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.mapcolfreq!-Tuple{Function,AbstractArray{MIToS.MSA.Residue,2},Union{MIToS.Information.Counts{T,1,A},MIToS.Information.Probabilities{T,1,A}}}",
    "page": "Information",
    "title": "MIToS.Information.mapcolfreq!",
    "category": "Method",
    "text": "It efficiently map a function (first argument) that takes a table of Counts or Probabilities (third argument). The table is filled in place with the counts or probabilities of each column from the msa (second argument).\n\nweights (default: NoClustering()): Weights to be used for table counting.\npseudocounts (default: NoPseudocount()): Pseudocount object to be applied to table.\npseudofrequencies (default: NoPseudofrequencies()): Pseudofrequencies to be applied to the normalized (probabilities) table.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.mapcolpairfreq!",
    "page": "Information",
    "title": "MIToS.Information.mapcolpairfreq!",
    "category": "Function",
    "text": "It efficiently map a function (first argument) that takes a table of Counts or Probabilities (third argument). The table is filled in place with the counts or probabilities of each pair of columns from the msa (second argument). The fourth positional argument usediagonal indicates if the function should be applied to identical element pairs (default to Val{true}).\n\nweights (default: NoClustering()): Weights to be used for table counting.\npseudocounts (default: NoPseudocount()): Pseudocount object to be applied to table.\npseudofrequencies (default: NoPseudofrequencies()): Pseudofrequencies to be applied to the normalized (probabilities) table.\ndiagonalvalue (default: 0): Value to fill diagonal elements if usediagonal is Val{false}.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.mapseqfreq!-Tuple{Function,AbstractArray{MIToS.MSA.Residue,2},Union{MIToS.Information.Counts{T,1,A},MIToS.Information.Probabilities{T,1,A}}}",
    "page": "Information",
    "title": "MIToS.Information.mapseqfreq!",
    "category": "Method",
    "text": "It efficiently map a function (first argument) that takes a table of Counts or Probabilities (third argument). The table is filled in place with the counts or probabilities of each sequence from the msa (second argument).\n\nweights (default: NoClustering()): Weights to be used for table counting.\npseudocounts (default: NoPseudocount()): Pseudocount object to be applied to table.\npseudofrequencies (default: NoPseudofrequencies()): Pseudofrequencies to be applied to the normalized (probabilities) table.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.mapseqpairfreq!",
    "page": "Information",
    "title": "MIToS.Information.mapseqpairfreq!",
    "category": "Function",
    "text": "It efficiently map a function (first argument) that takes a table of Counts or Probabilities (third argument). The table is filled in place with the counts or probabilities of each pair of sequences from the msa (second argument). The fourth positional argument usediagonal indicates if the function should be applied to identical element pairs (default to Val{true}).\n\nweights (default: NoClustering()): Weights to be used for table counting.\npseudocounts (default: NoPseudocount()): Pseudocount object to be applied to table.\npseudofrequencies (default: NoPseudofrequencies()): Pseudofrequencies to be applied to the normalized (probabilities) table.\ndiagonalvalue (default: 0): Value to fill diagonal elements if usediagonal is Val{false}.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.marginal_entropy-Tuple{Union{MIToS.Information.Counts{T,N,A},MIToS.Information.Probabilities{T,N,A}},Int64,Real}",
    "page": "Information",
    "title": "MIToS.Information.marginal_entropy",
    "category": "Method",
    "text": "It calculates marginal entropy (H) from a table of Counts or Probabilities. The second positional argument is used to indicate the magin used to calculate the entropy, e.g. it estimates the entropy H(X) if marginal is 1, H(Y) for 2, etc. Use last and optional positional argument to change the base of the log. The default base is e, so the result is in nats. You can use 2.0 as base to get the result in bits.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.mutual_information-Tuple{Union{MIToS.Information.Counts{T,N,A},MIToS.Information.Probabilities{T,N,A}},Real}",
    "page": "Information",
    "title": "MIToS.Information.mutual_information",
    "category": "Method",
    "text": "It calculates Mutual Information (MI) from a table of Counts or Probabilities. Use last and optional positional argument to change the base of the log. The default base is e, so the result is in nats. You can use 2.0 as base to get the result in bits. Calculation of MI from Counts is faster than from Probabilities.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.normalized_mutual_information-Tuple{Union{MIToS.Information.Counts{T,N,A},MIToS.Information.Probabilities{T,N,A}}}",
    "page": "Information",
    "title": "MIToS.Information.normalized_mutual_information",
    "category": "Method",
    "text": "It calculates a Normalized Mutual Information (nMI) by Entropy from a table of Counts or Probabilities.\n\nnMI(X, Y) = MI(X, Y) / H(X, Y)\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.pairwisegapfraction-Tuple{AbstractArray{MIToS.MSA.Residue,2}}",
    "page": "Information",
    "title": "MIToS.Information.pairwisegapfraction",
    "category": "Method",
    "text": "It takes a MSA or a file and a Format as first arguments. It calculates the percentage of gaps on columns pairs (union and intersection) using sequence clustering (Hobohm I).\n\nArgument, type, default value and descriptions:\n\n    - clustering  Bool      true    Sequence clustering (Hobohm I)\n    - threshold             62      Percent identity threshold for sequence clustering (Hobohm I)\n\nThis function returns:\n\n    - pairwise gap union as percentage\n    - pairwise gap intersection as percentage\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.probabilities!-Tuple{MIToS.Information.ContingencyTable{T,N,A},Any,MIToS.Information.Pseudocount,MIToS.Information.Pseudofrequencies,Vararg{AbstractArray{MIToS.MSA.Residue,1},N}}",
    "page": "Information",
    "title": "MIToS.Information.probabilities!",
    "category": "Method",
    "text": "It populates a ContingencyTable (first argument) using the probabilities in the sequences (last positional arguments). The dimension of the table must match the number of sequences and all the sequences must have the same length. You must indicate the used weights, pseudocounts and pseudofrequencies as second, third and fourth positional arguments respectively. You can use NoClustering(), NoPseudocount() and NoPseudofrequencies() to avoid the use of sequence weighting, pseudocounts and pseudofrequencies, respectively.\n\n\n\n"
},

{
    "location": "Information_API.html#MIToS.Information.probabilities-Tuple{Vararg{AbstractArray{MIToS.MSA.Residue,1},N}}",
    "page": "Information",
    "title": "MIToS.Information.probabilities",
    "category": "Method",
    "text": "It returns a ContingencyTable wrapped in a Probabilities type with the frequencies of residues in the sequences that takes as arguments. The dimension of the table is equal to the number of sequences. You can use the keyword arguments alphabet, weights, pseudocounts and pseudofrequencies to indicate the alphabet of the table (default to UngappedAlphabet()), a clustering result (default to NoClustering()),  the pseudocounts (default to NoPseudocount()) and the pseudofrequencies (default to NoPseudofrequencies()) to be used during the estimation of the probabilities.\n\n\n\n"
},

{
    "location": "Information_API.html#StatsBase.entropy-Tuple{Union{MIToS.Information.Counts{T,N,A},MIToS.Information.Probabilities{T,N,A}},Real}",
    "page": "Information",
    "title": "StatsBase.entropy",
    "category": "Method",
    "text": "It calculates the Shannon entropy (H) from a table of Counts or Probabilities. Use last and optional positional argument to change the base of the log. The default base is e, so the result is in nats. You can use 2.0 as base to get the result in bits.\n\n\n\n"
},

{
    "location": "Information_API.html#Methods-and-functions-1",
    "page": "Information",
    "title": "Methods and functions",
    "category": "section",
    "text": "Modules = [MIToS.Information]\nPrivate = false\nOrder   = [:function]"
},

{
    "location": "SIFTS_API.html#",
    "page": "SIFTS",
    "title": "SIFTS",
    "category": "page",
    "text": ""
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS",
    "page": "SIFTS",
    "title": "MIToS.SIFTS",
    "category": "Module",
    "text": "The SIFTS module of MIToS allows to obtain the residue-level mapping between databases stored in the SIFTS XML files. It makes easy to assign PDB residues to UniProt/Pfam positions. Given the fact that pairwise alignments can lead to misleading association between residues in both sequences, SIFTS offers  more reliable association between sequence and structure residue numbers.\n\nFeatures\n\nDownload and parse SIFTS XML files\nStore residue-level mapping in Julia\nEasy generation of Dicts between residues numbers\n\nusing MIToS.SIFTS\n\n\n\n"
},

{
    "location": "SIFTS_API.html#SIFTS-1",
    "page": "SIFTS",
    "title": "SIFTS",
    "category": "section",
    "text": "MIToS.SIFTS"
},

{
    "location": "SIFTS_API.html#Contents-1",
    "page": "SIFTS",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"SIFTS_API.md\"]\nDepth = 2"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.SIFTSResidue",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.SIFTSResidue",
    "category": "Type",
    "text": "A SIFTSResidue object stores the SIFTS residue level mapping for a residue. It has the following fields that you can access at any moment for query purposes:\n\n- `PDBe` : A `dbPDBe` object, it's present in all the `SIFTSResidue`s.\n- `UniProt` : A `Nullable` wrapping a `dbUniProt` object.\n- `Pfam` : A `Nullable` wrapping a `dbPfam` object.\n- `NCBI` : A `Nullable` wrapping a `dbNCBI` object.\n- `InterPro` : A `Nullable` wrapping a `dbInterPro` object.\n- `PDB` : A `Nullable` wrapping a `dbPDB` object.\n- `SCOP` : A `Nullable` wrapping a `dbSCOP` object.\n- `CATH` : A `Nullable` wrapping a `dbCATH` object.\n- `missing` : It's `true` if the residue is missing, i.e. not observed, in the structure.\n- `sscode` : A string with the secondary structure code of the residue.\n- `ssname` : A string with the secondary structure name of the residue.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.dbCATH",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.dbCATH",
    "category": "Type",
    "text": "dbCATH stores the residue id, number, name and chain in CATH as strings.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.dbInterPro",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.dbInterPro",
    "category": "Type",
    "text": "dbInterPro stores the residue id, number, name and evidence in InterPro as strings.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.dbNCBI",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.dbNCBI",
    "category": "Type",
    "text": "dbNCBI stores the residue id, number and name in NCBI as strings.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.dbPDB",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.dbPDB",
    "category": "Type",
    "text": "dbPDB stores the residue id, number, name and chain in PDB as strings.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.dbPDBe",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.dbPDBe",
    "category": "Type",
    "text": "dbPDBe stores the residue number and name in PDBe as strings.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.dbPfam",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.dbPfam",
    "category": "Type",
    "text": "dbPfam stores the residue id, number and name in Pfam as strings.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.dbSCOP",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.dbSCOP",
    "category": "Type",
    "text": "dbSCOP stores the residue id, number, name and chain in SCOP as strings.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.dbUniProt",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.dbUniProt",
    "category": "Type",
    "text": "dbUniProt stores the residue id, number and name in UniProt as strings.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#Types-1",
    "page": "SIFTS",
    "title": "Types",
    "category": "section",
    "text": "Modules = [MIToS.SIFTS]\nPrivate = false\nOrder   = [:type]"
},

{
    "location": "SIFTS_API.html#Constants-1",
    "page": "SIFTS",
    "title": "Constants",
    "category": "section",
    "text": "Modules = [MIToS.SIFTS]\nPrivate = false\nOrder   = [:constant]"
},

{
    "location": "SIFTS_API.html#Macros-1",
    "page": "SIFTS",
    "title": "Macros",
    "category": "section",
    "text": "Modules = [MIToS.SIFTS]\nPrivate = false\nOrder   = [:macro]"
},

{
    "location": "SIFTS_API.html#Base.parse-Tuple{LightXML.XMLDocument,Type{MIToS.SIFTS.SIFTSXML}}",
    "page": "SIFTS",
    "title": "Base.parse",
    "category": "Method",
    "text": "parse(document::LightXML.XMLDocument, ::Type{SIFTSXML}; chain=All, missings::Bool=true)\n\nReturns a Vector{SIFTSResidue} parsed from a SIFTSXML file. By default, parses all the chains and includes missings residues.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.downloadsifts-Tuple{String}",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.downloadsifts",
    "category": "Method",
    "text": "Download the gzipped SIFTS xml  for the pdbcode. The extension of the downloaded file is .xml.gz by default. The filename can be changed, but the .xml.gz at the end is mandatory.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#MIToS.SIFTS.siftsmapping-Tuple{String,Type{F},String,Type{T},String}",
    "page": "SIFTS",
    "title": "MIToS.SIFTS.siftsmapping",
    "category": "Method",
    "text": "Parses a SIFTS XML file and returns a Dict between residue numbers of two DataBases with the given identifiers. A chain could be specified (All by default). If missings is true (default) all the residues are used, even if they haven’t coordinates in the PDB file.\n\n\n\n"
},

{
    "location": "SIFTS_API.html#Methods-and-functions-1",
    "page": "SIFTS",
    "title": "Methods and functions",
    "category": "section",
    "text": "Modules = [MIToS.SIFTS]\nPrivate = false\nOrder   = [:function]"
},

{
    "location": "PDB_API.html#",
    "page": "PDB",
    "title": "PDB",
    "category": "page",
    "text": ""
},

{
    "location": "PDB_API.html#MIToS.PDB",
    "page": "PDB",
    "title": "MIToS.PDB",
    "category": "Module",
    "text": "The module PDB defines types and methods to work with protein structures inside Julia. It is useful to link structural and sequential information, and needed for measure the predictive performance at protein contact prediction of mutual information scores.\n\nFeatures\n\nRead and parse PDF and PDBML files\nCalculate distance and contacts between atoms or residues\nDetermine interaction between residues\n\nusing MIToS.PDB\n\n\n\n"
},

{
    "location": "PDB_API.html#PDB-1",
    "page": "PDB",
    "title": "PDB",
    "category": "section",
    "text": "MIToS.PDB"
},

{
    "location": "PDB_API.html#Contents-1",
    "page": "PDB",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"PDB_API.md\"]\nDepth = 2"
},

{
    "location": "PDB_API.html#MIToS.PDB.Coordinates",
    "page": "PDB",
    "title": "MIToS.PDB.Coordinates",
    "category": "Type",
    "text": "A Coordinates object is a fixed size vector with the coordinates x,y,z.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.PDBAtom",
    "page": "PDB",
    "title": "MIToS.PDB.PDBAtom",
    "category": "Type",
    "text": "A PDBAtom object contains the information from a PDB atom, without information of the residue. It has the following fields that you can access at any moment for query purposes:\n\n- `coordinates` : x,y,z coordinates, e.g. `Coordinates(109.641,73.162,42.7)`.\n- `atom` : Atom name, e.g. `\"CA\"`.\n- `element` : Element type of the atom, e.g. `\"C\"`.\n- `occupancy` : A float number with the occupancy, e.g. `1.0`.\n- `B` : B factor as a string, e.g. `\"23.60\"`.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.PDBFile",
    "page": "PDB",
    "title": "MIToS.PDB.PDBFile",
    "category": "Type",
    "text": "PDBFile <: Format\n\nProtein Data Bank (PDB) format. It provides a standard representation for macromolecular structure data derived from X-ray diffraction and NMR studies.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.PDBML",
    "page": "PDB",
    "title": "MIToS.PDB.PDBML",
    "category": "Type",
    "text": "PDBML <: Format\n\nProtein Data Bank Markup Language (PDBML), a representation of PDB data in XML format.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.PDBResidue",
    "page": "PDB",
    "title": "MIToS.PDB.PDBResidue",
    "category": "Type",
    "text": "A PDBResidue object contains all the information about a PDB residue. It has the following fields that you can access at any moment for query purposes:\n\n- `id` : A `PDBResidueIdentifier` object.\n- `atoms` : A vector of `PDBAtom`s.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.PDBResidueIdentifier",
    "page": "PDB",
    "title": "MIToS.PDB.PDBResidueIdentifier",
    "category": "Type",
    "text": "A PDBResidueIdentifier object contains the information needed to identity PDB residues. It has the following fields that you can access at any moment for query purposes:\n\n- `PDBe_number` : It's only used when a PDBML is readed (PDBe number as a string).\n- `number` : PDB residue number, it includes insertion codes, e.g. `\"34A\"`.\n- `name` : Three letter residue name in PDB, e.g. `\"LYS\"`.\n- `group` : It can be `\"ATOM\"` or `\"HETATM\"`.\n- `model` : The model number as a string, e.g. `\"1\"`.\n- `chain` : The chain as a string, e.g. `\"A\"`.\n\n\n\n"
},

{
    "location": "PDB_API.html#Types-1",
    "page": "PDB",
    "title": "Types",
    "category": "section",
    "text": "Modules = [MIToS.PDB]\nPrivate = false\nOrder   = [:type]"
},

{
    "location": "PDB_API.html#MIToS.PDB.covalentradius",
    "page": "PDB",
    "title": "MIToS.PDB.covalentradius",
    "category": "Constant",
    "text": "Covalent radius in Å of each element from the Additional file 1 of PICCOLO [1]. Hydrogen was updated using the value on Table 2 from Cordero et. al. [2].\n\nBickerton, G. R., Higueruelo, A. P., & Blundell, T. L. (2011).\n\nComprehensive, atomic-level characterization of structurally characterized protein-protein interactions: the PICCOLO database. BMC bioinformatics, 12(1), 313.\n\nCordero, B., Gómez, V., Platero-Prats, A. E., Revés, M.,\n\nEcheverría, J., Cremades, E., ... & Alvarez, S. (2008). Covalent radii revisited. Dalton Transactions, (21), 2832-2838.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.vanderwaalsradius",
    "page": "PDB",
    "title": "MIToS.PDB.vanderwaalsradius",
    "category": "Constant",
    "text": "van der Waals radius in Å from the Additional file 1 of Bickerton et. al. 2011\n\nBickerton, G. R., Higueruelo, A. P., & Blundell, T. L. (2011).\n\nComprehensive, atomic-level characterization of structurally characterized protein-protein interactions: the PICCOLO database. BMC bioinformatics, 12(1), 313.\n\n\n\n"
},

{
    "location": "PDB_API.html#Constants-1",
    "page": "PDB",
    "title": "Constants",
    "category": "section",
    "text": "Modules = [MIToS.PDB]\nPrivate = false\nOrder   = [:constant]"
},

{
    "location": "PDB_API.html#MIToS.PDB.@atoms-Tuple{Any,Symbol,Any,Symbol,Any,Symbol,Any,Symbol,Any,Symbol,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.@atoms",
    "category": "Macro",
    "text": "@atoms ... model ... chain ... group ... residue ... atom ...\n\nThese return a vector of PDBAtoms with the selected subset of atoms. You can use the type All to avoid filtering that option.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.@residues-Tuple{Any,Symbol,Any,Symbol,Any,Symbol,Any,Symbol,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.@residues",
    "category": "Macro",
    "text": "@residues ... model ... chain ... group ... residue ...\n\nThese return a new vector with the selected subset of residues from a list of residues. You can use the type All to avoid filtering that option.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.@residuesdict-Tuple{Any,Symbol,Any,Symbol,Any,Symbol,Any,Symbol,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.@residuesdict",
    "category": "Macro",
    "text": "@residuesdict ... model ... chain ... group ... residue ...\n\nThese return a dictionary (using PDB residue numbers as keys) with the selected subset of residues. You can use the type All to avoid filtering that option.\n\n\n\n"
},

{
    "location": "PDB_API.html#Macros-1",
    "page": "PDB",
    "title": "Macros",
    "category": "section",
    "text": "Modules = [MIToS.PDB]\nPrivate = false\nOrder   = [:macro]"
},

{
    "location": "PDB_API.html#Base.angle-Tuple{MIToS.PDB.Coordinates,MIToS.PDB.Coordinates,MIToS.PDB.Coordinates}",
    "page": "PDB",
    "title": "Base.angle",
    "category": "Method",
    "text": "angle(a::Coordinates, b::Coordinates, c::Coordinates)\n\nAngle (in degrees) at b between a-b and b-c\n\n\n\n"
},

{
    "location": "PDB_API.html#Base.any-Tuple{Function,MIToS.PDB.PDBResidue,MIToS.PDB.PDBResidue,Function}",
    "page": "PDB",
    "title": "Base.any",
    "category": "Method",
    "text": "any(f::Function, a::PDBResidue, b::PDBResidue, criteria::Function)\n\nTest if the function f is true for any pair of atoms between the residues a and b. This function only test atoms that returns true for the fuction criteria.\n\n\n\n"
},

{
    "location": "PDB_API.html#Base.any-Tuple{Function,MIToS.PDB.PDBResidue,MIToS.PDB.PDBResidue}",
    "page": "PDB",
    "title": "Base.any",
    "category": "Method",
    "text": "any(f::Function, a::PDBResidue, b::PDBResidue)\n\nTest if the function f is true for any pair of atoms between the residues a and b\n\n\n\n"
},

{
    "location": "PDB_API.html#Base.parse-Tuple{LightXML.XMLDocument,Type{MIToS.PDB.PDBML}}",
    "page": "PDB",
    "title": "Base.parse",
    "category": "Method",
    "text": "parse(pdbml, ::Type{PDBML}; chain=All, model=All, group=All, atomname=All, onlyheavy=false, label=true, occupancyfilter=false)\n\nReads a LightXML.XMLDocument representing a pdb file. Returns a list of PDBResidues (view MIToS.PDB.PDBResidues). Setting chain, model, group, atomname and onlyheavy values can be used to select of a subset of all residues. If not set, all residues are returned. If the keyword argument label (default: true) is false,the auth_ attributes will be use instead of the label_ attributes for chain, atom and residue name fields. The auth_ attributes are alternatives provided by an author in order to match the identification/values used in the publication that describes the structure. If the keyword argument occupancyfilter (default: false) is true, only the atoms with the best occupancy are returned.\n\n\n\n"
},

{
    "location": "PDB_API.html#Base.parse-Tuple{Union{IO,String},Type{MIToS.PDB.PDBFile}}",
    "page": "PDB",
    "title": "Base.parse",
    "category": "Method",
    "text": "parse(io, ::Type{PDBFile}; chain=All, model=All, group=All, atomname=All, onlyheavy=false, occupancyfilter=false)\n\nReads a text file of a PDB entry. Returns a list of PDBResidue (view MIToS.PDB.PDBResidues). Setting chain, model, group, atomname and onlyheavy values can be used to select of a subset of all residues. Group can be \"ATOM\" or \"HETATM\". If not set, all residues are returned. If the keyword argument occupancyfilter (default: false) is true, only the atoms with the best occupancy are returned.\n\n\n\n"
},

{
    "location": "PDB_API.html#Base.print",
    "page": "PDB",
    "title": "Base.print",
    "category": "Function",
    "text": "print(io, res, format::Type{PDBFile}) print(res, format::Type{PDBFile})\n\nPrint a PDBResidue or a vector of PDBResidues in PDB format.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.CAmatrix-Tuple{AbstractArray{MIToS.PDB.PDBResidue,1}}",
    "page": "PDB",
    "title": "MIToS.PDB.CAmatrix",
    "category": "Method",
    "text": "Returns a matrix with the x, y and z coordinates of the Cα with best occupancy for each PDBResidue.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.aromatic-Tuple{MIToS.PDB.PDBResidue,MIToS.PDB.PDBResidue}",
    "page": "PDB",
    "title": "MIToS.PDB.aromatic",
    "category": "Method",
    "text": "There's an aromatic interaction if centriods are at 6.0 Å or less.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.aromaticsulphur-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.PDBAtom,Any,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.aromaticsulphur",
    "category": "Method",
    "text": "Returns true if an sulphur and an aromatic atoms are 5.3 Å or less\"\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.atoms",
    "page": "PDB",
    "title": "MIToS.PDB.atoms",
    "category": "Function",
    "text": "atoms(residue_list, model, chain, group, residue, atom)\n\nThese return a vector of PDBAtoms with the selected subset of atoms. You can use the type All (default value of the positional arguments) to avoid filtering a that level.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.bestoccupancy-Tuple{Array{MIToS.PDB.PDBAtom,1}}",
    "page": "PDB",
    "title": "MIToS.PDB.bestoccupancy",
    "category": "Method",
    "text": "Takes a Vector of PDBAtoms and returns a Vector of the PDBAtoms with best occupancy.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.center!-Tuple{Array{Float64,2}}",
    "page": "PDB",
    "title": "MIToS.PDB.center!",
    "category": "Method",
    "text": "center!(A::Matrix{Float64})\n\nTakes a set of points A as an NxD matrix (N: number of points, D: dimension). Translates A in place so that its centroid is at the origin of coordinates\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.centeredcoordinates",
    "page": "PDB",
    "title": "MIToS.PDB.centeredcoordinates",
    "category": "Function",
    "text": "Returns a Matrix{Float64} with the centered coordinates of all the atoms in residues. An optional positional argument CA (default: true) defines if only Cα carbons should be used to center the matrix.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.centeredresidues",
    "page": "PDB",
    "title": "MIToS.PDB.centeredresidues",
    "category": "Function",
    "text": "Returns a new Vector{PDBResidue} with the PDBResidues having centered coordinates. An optional positional argument CA (default: true) defines if only Cα carbons should be used to center the matrix.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.change_coordinates",
    "page": "PDB",
    "title": "MIToS.PDB.change_coordinates",
    "category": "Function",
    "text": "change_coordinates(residue::PDBResidue, coordinates::Matrix{Float64}, offset::Int=1)\n\nReturns a new PDBResidues with (x,y,z) from a coordinates Matrix{Float64} You can give an offset indicating in wich matrix row starts the (x,y,z) coordinates of the residue.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.change_coordinates-Tuple{AbstractArray{MIToS.PDB.PDBResidue,1},Array{Float64,2}}",
    "page": "PDB",
    "title": "MIToS.PDB.change_coordinates",
    "category": "Method",
    "text": "change_coordinates(residues::AbstractVector{PDBResidue}, coordinates::Matrix{Float64})\n\nReturns a new Vector{PDBResidues} with (x,y,z) from a coordinates Matrix{Float64}\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.change_coordinates-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.Coordinates}",
    "page": "PDB",
    "title": "MIToS.PDB.change_coordinates",
    "category": "Method",
    "text": "change_coordinates(atom::PDBAtom, coordinates::Coordinates)\n\nReturns a new PDBAtom but with a new coordinates\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.check_atoms_for_interactions-Tuple{MIToS.PDB.PDBResidue}",
    "page": "PDB",
    "title": "MIToS.PDB.check_atoms_for_interactions",
    "category": "Method",
    "text": "This function takes a PDBResidue and returns true only if all the atoms can be used for checking interactions.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.contact-Tuple{Array{MIToS.PDB.PDBResidue,1},AbstractFloat}",
    "page": "PDB",
    "title": "MIToS.PDB.contact",
    "category": "Method",
    "text": "contact(residues::Vector{PDBResidue}, limit::AbstractFloat; criteria::String=\"All\")\n\nIf contact takes a Vector{PDBResidue}, It returns a matrix with all the pairwise comparisons (contact map).\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.contact-Tuple{MIToS.PDB.Coordinates,MIToS.PDB.Coordinates,AbstractFloat}",
    "page": "PDB",
    "title": "MIToS.PDB.contact",
    "category": "Method",
    "text": "contact(a::Coordinates, b::Coordinates, limit::AbstractFloat)\n\nIt returns true if the distance is less or equal to the limit. It doesn't call sqrt because it does squared_distance(a,b) <= limit^2.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.contact-Tuple{MIToS.PDB.PDBResidue,MIToS.PDB.PDBResidue,AbstractFloat}",
    "page": "PDB",
    "title": "MIToS.PDB.contact",
    "category": "Method",
    "text": "contact(A::PDBResidue, B::PDBResidue, limit::AbstractFloat; criteria::String=\"All\")\n\nReturns true if the residues A and B are at contact distance (limit). The available distance criteria are: Heavy, All, CA, CB (CA for GLY)\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.coordinatesmatrix-Tuple{MIToS.PDB.PDBResidue}",
    "page": "PDB",
    "title": "MIToS.PDB.coordinatesmatrix",
    "category": "Method",
    "text": "Returns a matrix with the x, y, z coordinates of each atom in each PDBResidue\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.covalent-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.PDBAtom,Any,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.covalent",
    "category": "Method",
    "text": "Returns true if the distance between atoms is less than the sum of the covalentradius of each atom.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.distance-Tuple{Array{MIToS.PDB.PDBResidue,1}}",
    "page": "PDB",
    "title": "MIToS.PDB.distance",
    "category": "Method",
    "text": "distance(residues::Vector{PDBResidue}; criteria::String=\"All\")\n\nIf distance takes a Vector{PDBResidue} returns a PairwiseListMatrix{Float64, false} with all the pairwise comparisons (distance matrix).\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.distance-Tuple{MIToS.PDB.Coordinates,MIToS.PDB.Coordinates}",
    "page": "PDB",
    "title": "MIToS.PDB.distance",
    "category": "Method",
    "text": "It calculates the squared euclidean distance.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.disulphide-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.PDBAtom,Any,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.disulphide",
    "category": "Method",
    "text": "Returns true if two CYS's S are at 2.08 Å or less\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.downloadpdb-Tuple{String}",
    "page": "PDB",
    "title": "MIToS.PDB.downloadpdb",
    "category": "Method",
    "text": "Download a gzipped PDB file from PDB database. Requires a four character pdbcode. By default the format is PDBML (PDB XML) and uses the baseurl http://www.rcsb.org/pdb/files/. filename is the path/name of the output file.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.findCB-Tuple{MIToS.PDB.PDBResidue}",
    "page": "PDB",
    "title": "MIToS.PDB.findCB",
    "category": "Method",
    "text": "Returns a vector of indices for CB (CA for GLY)\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.findatoms-Tuple{Array{MIToS.PDB.PDBAtom,1},String}",
    "page": "PDB",
    "title": "MIToS.PDB.findatoms",
    "category": "Method",
    "text": "findatoms(res::PDBResidue, atom::String)\n\nReturns a index vector of the atoms with the given atom name.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.findheavy-Tuple{Array{MIToS.PDB.PDBAtom,1}}",
    "page": "PDB",
    "title": "MIToS.PDB.findheavy",
    "category": "Method",
    "text": "Returns a list with the index of the heavy atoms (all atoms except hydrogen) in the PDBResidue\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.getCA-Tuple{MIToS.PDB.PDBResidue}",
    "page": "PDB",
    "title": "MIToS.PDB.getCA",
    "category": "Method",
    "text": "Returns the Cα with best occupancy in the PDBResidue.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.getpdbdescription-Tuple{String}",
    "page": "PDB",
    "title": "MIToS.PDB.getpdbdescription",
    "category": "Method",
    "text": "Access general information about a PDB entry (e.g., Header information) using the RESTful interface of the PDB database (describePDB). Returns a Dict for the four character pdbcode.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.hydrogenbond-Tuple{MIToS.PDB.PDBResidue,MIToS.PDB.PDBResidue}",
    "page": "PDB",
    "title": "MIToS.PDB.hydrogenbond",
    "category": "Method",
    "text": "This function only works if there are hydrogens in the structure. The criteria for a hydrogen bond are:\n\nd(Ai, Aj) < 3.9Å\nd(Ah, Aacc) < 2.5Å\nθ(Adon, Ah, Aacc) > 90°\nθ(Adon, Aacc, Aacc-antecedent) > 90°\nθ(Ah, Aacc, Aacc-antecedent) > 90°\n\nWhere Ah is the donated hydrogen atom, Adon is the hydrogen bond donor atom, Aacc is the hydrogen bond acceptor atom and Aacc-antecednt is the atom antecedent to the hydrogen bond acceptor atom.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.hydrophobic-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.PDBAtom,Any,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.hydrophobic",
    "category": "Method",
    "text": "There's an hydrophobic interaction if two hydrophobic atoms are at 5.0 Å or less.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.ionic-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.PDBAtom,Any,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.ionic",
    "category": "Method",
    "text": "There's an ionic interaction if a cationic and an anionic atoms are at 6.0 Å or less.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.isanionic-Tuple{MIToS.PDB.PDBAtom,String}",
    "page": "PDB",
    "title": "MIToS.PDB.isanionic",
    "category": "Method",
    "text": "Returns true if the atom, e.g. (\"GLU\",\"CD\"), is an anionic atom in the residue.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.isaromatic-Tuple{MIToS.PDB.PDBAtom,String}",
    "page": "PDB",
    "title": "MIToS.PDB.isaromatic",
    "category": "Method",
    "text": "Returns true if the atom, e.g. (\"HIS\",\"CG\"), is an aromatic atom in the residue.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.isatom-Tuple{MIToS.PDB.PDBAtom,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.isatom",
    "category": "Method",
    "text": "It tests if the atom has the indicated atom name.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.iscationic-Tuple{MIToS.PDB.PDBAtom,String}",
    "page": "PDB",
    "title": "MIToS.PDB.iscationic",
    "category": "Method",
    "text": "Returns true if the atom, e.g. (\"ARG\",\"NE\"), is a cationic atom in the residue.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.ishbondacceptor-Tuple{MIToS.PDB.PDBAtom,String}",
    "page": "PDB",
    "title": "MIToS.PDB.ishbondacceptor",
    "category": "Method",
    "text": "Returns true if the atom, e.g. (\"ARG\",\"O\"), is an acceptor in H bonds.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.ishbonddonor-Tuple{MIToS.PDB.PDBAtom,String}",
    "page": "PDB",
    "title": "MIToS.PDB.ishbonddonor",
    "category": "Method",
    "text": "Returns true if the atom, e.g. (\"ARG\",\"N\"), is a donor in H bonds.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.isresidue-Tuple{MIToS.PDB.PDBResidueIdentifier,Any,Any,Any,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.isresidue",
    "category": "Method",
    "text": "It tests if the PDB residue has the indicated model, chain, group and residue number.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.kabsch-Tuple{Array{Float64,2},Array{Float64,2}}",
    "page": "PDB",
    "title": "MIToS.PDB.kabsch",
    "category": "Method",
    "text": "kabsch(A::Matrix{Float64}, B::Matrix{Float64})\n\nThis function takes two sets of points, A (refrence) and B as NxD matrices, where D is the dimension and N is the number of points. Assumes that the centroids of A and B are at the origin of coordinates. You can call center! on each matrix before calling kabsch to center the matrices in the (0.0, 0.0, 0.0). Rotates B so that rmsd(A,B) is minimized. Returns the rotation matrix. You should do B * RotationMatrix to get the rotated B.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.mean_coordinates-Tuple{AbstractArray{Array{Float64,2},1}}",
    "page": "PDB",
    "title": "MIToS.PDB.mean_coordinates",
    "category": "Method",
    "text": "Calculates the average/mean position of each atom in a set of structure. The function takes a vector (AbstractVector) of vectors (Vector{PDBResidue}) or matrices (Matrix{Float64}) as first argument. As second (optional) argument this function can take an AbstractVector{Float64} of matrix/structure weights to return a weighted mean. When a Vector{PDBResidue} is used, if the keyword argument calpha is false the RMSF is calculated for all the atoms. By default only alpha carbons are used (default: calpha=true).\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.pication-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.PDBAtom,Any,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.pication",
    "category": "Method",
    "text": "There's a Π-Cation interaction if a cationic and an aromatic atoms are at 6.0 Å or less\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.proximitymean",
    "page": "PDB",
    "title": "MIToS.PDB.proximitymean",
    "category": "Function",
    "text": "proximitymean calculates the proximity mean/average for each residue as the average score (from a scores list) of all the residues within a certain physical distance to a given amino acid. The score of that residue is not included in the mean unless you set include to true. The default values are 6.05 for the distance threshold/limit and \"Heavy\" for the criteria keyword argument. This function allows to calculate pMI (proximity mutual information) and pC (proximity conservation) as in Buslje et. al. 2010.\n\nBuslje, Cristina Marino, Elin Teppa, Tomas Di Doménico, José María Delfino, and Morten Nielsen. Networks of high mutual information define the structural proximity of catalytic sites: implications for catalytic residue identification. PLoS Comput Biol 6, no. 11 (2010): e1000978.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.residuepairsmatrix-Tuple{Array{MIToS.PDB.PDBResidue,1},Type{T},Type{Val{diagonal}},T}",
    "page": "PDB",
    "title": "MIToS.PDB.residuepairsmatrix",
    "category": "Method",
    "text": "It creates a NamedArray containing a PairwiseListMatrix where each element (column, row) is identified with a PDBResidue from the input vector. You can indicate the value type of the matrix (default to Float64), if the list should have the diagonal values (default to Val{false}) and the diagonal values (default to NaN).\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.residues",
    "page": "PDB",
    "title": "MIToS.PDB.residues",
    "category": "Function",
    "text": "residues(residue_list, model, chain, group, residue)\n\nThese return a new vector with the selected subset of residues from a list of residues. You can use the type All (default value) to avoid filtering a that level.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.residuesdict",
    "page": "PDB",
    "title": "MIToS.PDB.residuesdict",
    "category": "Function",
    "text": "residuesdict(residue_list, model, chain, group, residue)\n\nThese return a dictionary (using PDB residue numbers as keys) with the selected subset of residues. You can use the type All (default value) to avoid filtering a that level.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.rmsd-Tuple{AbstractArray{MIToS.PDB.PDBResidue,1},AbstractArray{MIToS.PDB.PDBResidue,1}}",
    "page": "PDB",
    "title": "MIToS.PDB.rmsd",
    "category": "Method",
    "text": "rmsd(A::AbstractVector{PDBResidue}, B::AbstractVector{PDBResidue}; superimposed::Bool=false)\n\nReturns the Cα RMSD value between two PDB structures: A and B. If the structures are already superimposed between them, use superimposed=true to avoid a new superimposition (superimposed is false by default).\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.rmsd-Tuple{Array{Float64,2},Array{Float64,2}}",
    "page": "PDB",
    "title": "MIToS.PDB.rmsd",
    "category": "Method",
    "text": "rmsd(A::Matrix{Float64}, B::Matrix{Float64})\n\nReturn RMSD between two sets of points A and B, given as NxD matrices (N: number of points, D: dimension).\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.rmsf-Tuple{AbstractArray{Array{Float64,2},1}}",
    "page": "PDB",
    "title": "MIToS.PDB.rmsf",
    "category": "Method",
    "text": "Calculates the RMSF (Root Mean-Square-Fluctuation) between an atom and its average position in a set of structures. The function takes a vector (AbstractVector) of vectors (Vector{PDBResidue}) or matrices (Matrix{Float64}) as first argument. As second (optional) argument this function can take an AbstractVector{Float64} of matrix/structure weights to return the root weighted mean-square-fluctuation around the weighted mean structure. When a Vector{PDBResidue} is used, if the keyword argument calpha is false the RMSF is calculated for all the atoms. By default only alpha carbons are used (default: calpha=true).\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.selectbestoccupancy-Tuple{Array{MIToS.PDB.PDBAtom,1},Array{Int64,1}}",
    "page": "PDB",
    "title": "MIToS.PDB.selectbestoccupancy",
    "category": "Method",
    "text": "Takes a PDBResidue and a Vector of atom indices. Returns the index value of the Vector with maximum occupancy.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.squared_distance-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.PDBAtom}",
    "page": "PDB",
    "title": "MIToS.PDB.squared_distance",
    "category": "Method",
    "text": "It calculates the squared euclidean distance, i.e. it doesn't spend time in sqrt\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.squared_distance-Tuple{MIToS.PDB.PDBResidue,MIToS.PDB.PDBResidue}",
    "page": "PDB",
    "title": "MIToS.PDB.squared_distance",
    "category": "Method",
    "text": "squared_distance(A::PDBResidue, B::PDBResidue; criteria::String=\"All\")\n\nReturns the squared distance between the residues A and B. The available criteria are: Heavy, All, CA, CB (CA for GLY)\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.superimpose-Tuple{Array{MIToS.PDB.PDBResidue,1},Array{MIToS.PDB.PDBResidue,1}}",
    "page": "PDB",
    "title": "MIToS.PDB.superimpose",
    "category": "Method",
    "text": "This function takes A::Vector{PDBResidue} (reference) and B::Vector{PDBResidue}. Translates A and B to the origin of coordinates, and rotates B so that rmsd(A,B) is minimized with the Kabsch algorithm (using only their α carbons). Returns the rotated and translated versions of A and B, and the RMSD value.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.vanderwaals-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.PDBAtom,Any,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.vanderwaals",
    "category": "Method",
    "text": "Test if two atoms or residues are in van der Waals contact using: distance(a,b) <= 0.5 + vanderwaalsradius[a] + vanderwaalsradius[b]. It returns distance <= 0.5 if the atoms aren't in vanderwaalsradius.\n\n\n\n"
},

{
    "location": "PDB_API.html#MIToS.PDB.vanderwaalsclash-Tuple{MIToS.PDB.PDBAtom,MIToS.PDB.PDBAtom,Any,Any}",
    "page": "PDB",
    "title": "MIToS.PDB.vanderwaalsclash",
    "category": "Method",
    "text": "Returns true if the distance between the atoms is less than the sum of the vanderwaalsradius of the atoms. If the atoms aren't on the list (i.e. OXT), the vanderwaalsradius of the element is used. If there is not data in the dict, distance 0.0 is used.\n\n\n\n"
},

{
    "location": "PDB_API.html#Methods-and-functions-1",
    "page": "PDB",
    "title": "Methods and functions",
    "category": "section",
    "text": "Modules = [MIToS.PDB]\nPrivate = false\nOrder   = [:function]"
},

{
    "location": "Pfam_API.html#",
    "page": "Pfam",
    "title": "Pfam",
    "category": "page",
    "text": ""
},

{
    "location": "Pfam_API.html#MIToS.Pfam",
    "page": "Pfam",
    "title": "MIToS.Pfam",
    "category": "Module",
    "text": "The Pfam module, defines functions to measure the protein contact prediction performance of information measure between column pairs from a Pfam MSA.\n\nFeatures\n\nRead and download Pfam MSAs\nObtain PDB information from alignment annotations\nMap between sequence/alignment residues/columns and PDB structures\nMeasure of AUC (ROC curve) for contact prediction of MI scores\n\nusing MIToS.Pfam\n\n\n\n"
},

{
    "location": "Pfam_API.html#Pfam-1",
    "page": "Pfam",
    "title": "Pfam",
    "category": "section",
    "text": "MIToS.Pfam"
},

{
    "location": "Pfam_API.html#Contents-1",
    "page": "Pfam",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"Pfam_API.md\"]\nDepth = 2"
},

{
    "location": "Pfam_API.html#Types-1",
    "page": "Pfam",
    "title": "Types",
    "category": "section",
    "text": "Modules = [MIToS.Pfam]\nPrivate = false\nOrder   = [:type]"
},

{
    "location": "Pfam_API.html#Constants-1",
    "page": "Pfam",
    "title": "Constants",
    "category": "section",
    "text": "Modules = [MIToS.Pfam]\nPrivate = false\nOrder   = [:constant]"
},

{
    "location": "Pfam_API.html#Macros-1",
    "page": "Pfam",
    "title": "Macros",
    "category": "section",
    "text": "Modules = [MIToS.Pfam]\nPrivate = false\nOrder   = [:macro]"
},

{
    "location": "Pfam_API.html#MIToS.Pfam.downloadpfam-Tuple{String}",
    "page": "Pfam",
    "title": "MIToS.Pfam.downloadpfam",
    "category": "Method",
    "text": "Download a gzipped stockholm full alignment for the pfamcode. The extension of the downloaded file is .stockholm.gz by default. The filename can be changed, but the .gz at the end is mandatory.\n\n\n\n"
},

{
    "location": "Pfam_API.html#MIToS.Pfam.getcontactmasks-Tuple{Array{T<:AbstractFloat,1}}",
    "page": "Pfam",
    "title": "MIToS.Pfam.getcontactmasks",
    "category": "Method",
    "text": "This function takes a msacontacts or its list of contacts contact_list with 1.0 for true contacts and 0.0 for not contacts (NaN or other numbers for missing values). Returns two BitVectors, the first with trues where contact_list is 1.0 and the second with trues where contact_list is 0.0. There are useful for AUC calculations.\n\n\n\n"
},

{
    "location": "Pfam_API.html#MIToS.Pfam.getseq2pdb-Tuple{MIToS.MSA.AnnotatedMultipleSequenceAlignment}",
    "page": "Pfam",
    "title": "MIToS.Pfam.getseq2pdb",
    "category": "Method",
    "text": "Generates from a Pfam msa a Dict{String, Vector{Tuple{String,String}}}. Keys are sequence IDs and each value is a list of tuples containing PDB code and chain.\n\njulia> getseq2pdb(msa)\nDict{String,Array{Tuple{String,String},1}} with 1 entry:\n  \"F112_SSV1/3-112\" => [(\"2VQC\",\"A\")]\n\n\n\n\n"
},

{
    "location": "Pfam_API.html#MIToS.Pfam.hasresidues-Tuple{MIToS.MSA.AnnotatedMultipleSequenceAlignment,Dict{Int64,String}}",
    "page": "Pfam",
    "title": "MIToS.Pfam.hasresidues",
    "category": "Method",
    "text": "Returns a BitVector where there is a true for each column with PDB residue.\n\n\n\n"
},

{
    "location": "Pfam_API.html#MIToS.Pfam.msacolumn2pdbresidue-Tuple{MIToS.MSA.AnnotatedMultipleSequenceAlignment,String,String,String,String,String}",
    "page": "Pfam",
    "title": "MIToS.Pfam.msacolumn2pdbresidue",
    "category": "Method",
    "text": "msacolumn2pdbresidue(msa, seqid, pdbid, chain, pfamid, siftsfile; strict=false, checkpdbname=false, missings=true)\n\nThis function returns a Dict{Int,String} with MSA column numbers on the input file as keys and PDB residue numbers (\"\" for missings) as values. The mapping is performed using SIFTS. This function needs correct ColMap and SeqMap annotations. This checks correspondence of the residues between the MSA sequence and SIFTS (It throws a warning if there are differences). Missing residues are included if the keyword argument missings is true (default: true). If the keyword argument strict is true (default: false), throws an Error, instead of a Warning, when residues don't match. If the keyword argument checkpdbname is true (default: false), throws an Error if the three letter name of the PDB residue isn't the MSA residue. If you are working with a downloaded Pfam MSA without modifications, you should read it using generatemapping=true and useidcoordinates=true. If you don't indicate the path to the siftsfile used in the mapping, this function downloads the SIFTS file in the current folder. If you don't indicate the Pfam accession number (pfamid), this function tries to read the AC file annotation.\n\n\n\n"
},

{
    "location": "Pfam_API.html#MIToS.Pfam.msacontacts",
    "page": "Pfam",
    "title": "MIToS.Pfam.msacontacts",
    "category": "Function",
    "text": "This function takes an AnnotatedMultipleSequenceAlignment with correct ColMap annotations and two dicts:\n\nThe first is an OrderedDict{String,PDBResidue} from PDB residue number to PDBResidue.\nThe second is a Dict{Int,String} from MSA column number on the input file to PDB residue number.\n\nmsacontacts returns a PairwiseListMatrix{Float64,false} of 0.0 and 1.0 where 1.0 indicates a residue contact. Contacts are defined with an inter residue distance less or equal to distance_limit (default to 6.05) angstroms between any heavy atom. NaN indicates a missing value.\n\n\n\n"
},

{
    "location": "Pfam_API.html#MIToS.Pfam.msaresidues-Tuple{MIToS.MSA.AnnotatedMultipleSequenceAlignment,DataStructures.OrderedDict{String,MIToS.PDB.PDBResidue},Dict{Int64,String}}",
    "page": "Pfam",
    "title": "MIToS.Pfam.msaresidues",
    "category": "Method",
    "text": "This function takes an AnnotatedMultipleSequenceAlignment with correct ColMap annotations and two dicts:\n\nThe first is an OrderedDict{String,PDBResidue} from PDB residue number to PDBResidue.\nThe second is a Dict{Int,String} from MSA column number on the input file to PDB residue number.\n\nmsaresidues returns an OrderedDict{Int,PDBResidue} from input column number (ColMap) to PDBResidue. Residues on inserts are not included.\n\n\n\n"
},

{
    "location": "Pfam_API.html#ROCAnalysis.AUC-Tuple{Array{T,1},BitArray{1},BitArray{1}}",
    "page": "Pfam",
    "title": "ROCAnalysis.AUC",
    "category": "Method",
    "text": "AUC(scores_list::Vector, true_contacts::BitVector, false_contacts::BitVector)\n\nReturns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the scores_list for true_contacts prediction. The three vectors should have the same length and false_contacts should be true where there are not contacts.\n\n\n\n"
},

{
    "location": "Pfam_API.html#ROCAnalysis.AUC-Tuple{NamedArrays.NamedArray{L,2,PairwiseListMatrices.PairwiseListMatrix{L,false,VL},NL},BitArray{1},BitArray{1}}",
    "page": "Pfam",
    "title": "ROCAnalysis.AUC",
    "category": "Method",
    "text": "AUC(scores::PairwiseListMatrix, true_contacts::BitVector, false_contacts::BitVector)\n\nReturns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the scores for true_contacts prediction. scores, true_contacts and false_contacts should have the same number of elements and false_contacts should be true where there are not contacts.\n\n\n\n"
},

{
    "location": "Pfam_API.html#ROCAnalysis.AUC-Tuple{NamedArrays.NamedArray{L<:AbstractFloat,2,PairwiseListMatrices.PairwiseListMatrix{L<:AbstractFloat,false,VL},NL},NamedArrays.NamedArray{R<:AbstractFloat,2,PairwiseListMatrices.PairwiseListMatrix{L<:AbstractFloat,false,VR},NR}}",
    "page": "Pfam",
    "title": "ROCAnalysis.AUC",
    "category": "Method",
    "text": "AUC(scores::PairwiseListMatrix, msacontacts::PairwiseListMatrix)\n\nReturns the Area Under a ROC (Receiver Operating Characteristic) Curve (AUC) of the scores for msacontact prediction. score and msacontact lists are vinculated (inner join) by their labels (i.e. column number in the file). msacontact should have 1.0 for true contacts and 0.0 for not contacts (NaN or other numbers for missing values).\n\n\n\n"
},

{
    "location": "Pfam_API.html#Methods-and-functions-1",
    "page": "Pfam",
    "title": "Methods and functions",
    "category": "section",
    "text": "Modules = [MIToS.Pfam]\nPrivate = false\nOrder   = [:function]"
},

{
    "location": "Utils_API.html#",
    "page": "Utils",
    "title": "Utils",
    "category": "page",
    "text": ""
},

{
    "location": "Utils_API.html#MIToS.Utils",
    "page": "Utils",
    "title": "MIToS.Utils",
    "category": "Module",
    "text": "The Utils has common utils functions and types used in other modules.\n\nusing MIToS.Utils\n\n\n\n"
},

{
    "location": "Utils_API.html#API-Utils-1",
    "page": "Utils",
    "title": "Utils",
    "category": "section",
    "text": "MIToS.Utils"
},

{
    "location": "Utils_API.html#Contents-1",
    "page": "Utils",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"Utils_API.md\"]\nDepth = 2"
},

{
    "location": "Utils_API.html#MIToS.Utils.All",
    "page": "Utils",
    "title": "MIToS.Utils.All",
    "category": "Type",
    "text": "All is used instead of MIToS 1.0 \"all\" or \"*\", because it's possible to dispatch on it.\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.Format",
    "page": "Utils",
    "title": "MIToS.Utils.Format",
    "category": "Type",
    "text": "Format is used for write special parse (and read) methods on it.\n\n\n\n"
},

{
    "location": "Utils_API.html#Types-1",
    "page": "Utils",
    "title": "Types",
    "category": "section",
    "text": "Modules = [MIToS.Utils]\nPrivate = false\nOrder   = [:type]"
},

{
    "location": "Utils_API.html#Constants-1",
    "page": "Utils",
    "title": "Constants",
    "category": "section",
    "text": "Modules = [MIToS.Utils]\nPrivate = false\nOrder   = [:constant]"
},

{
    "location": "Utils_API.html#Macros-1",
    "page": "Utils",
    "title": "Macros",
    "category": "section",
    "text": "Modules = [MIToS.Utils]\nPrivate = false\nOrder   = [:macro]"
},

{
    "location": "Utils_API.html#Base.read-Tuple{AbstractString,Type{T<:MIToS.Utils.Format},Vararg{Any,N}}",
    "page": "Utils",
    "title": "Base.read",
    "category": "Method",
    "text": "read(pathname, Format [, Type [, … ] ] ) -> Type\n\nThis function opens a file in the pathname and calls parse(io, ...) for the given Format and Type on it. If the  pathname is an HTTP or FTP URL, the file is downloaded with download in a temporal file. Gzipped files should end on .gz.\n\n\n\n"
},

{
    "location": "Utils_API.html#Base.write",
    "page": "Utils",
    "title": "Base.write",
    "category": "Function",
    "text": "write{T<:Format}(filename::AbstractString, object, format::Type{T}, mode::ASCIIString=\"w\")\n\nThis function opens a file with filename and mode (default: \"w\") and writes (print) the object with the given format. Gzipped files should end on .gz.\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.check_file-Tuple{Any}",
    "page": "Utils",
    "title": "MIToS.Utils.check_file",
    "category": "Method",
    "text": "Returns the filename. Throws an ErrorException if the file doesn't exist, or a warning if the file is empty.\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.check_pdbcode-Tuple{String}",
    "page": "Utils",
    "title": "MIToS.Utils.check_pdbcode",
    "category": "Method",
    "text": "It checks if a PDB code has the correct format.\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.download_file-Tuple{AbstractString,AbstractString}",
    "page": "Utils",
    "title": "MIToS.Utils.download_file",
    "category": "Method",
    "text": "download_file uses Requests.jl instead of system calls to download files from the web. It takes the file url as first argument and, optionally, a path to save it. Keyword arguments (ie. allow_redirects, max_redirects, timeout, headers) are are directly passed to to Requests.get_streaming.\n\njulia> download_file(\"http://www.uniprot.org/uniprot/P69905.fasta\",\"seq.fasta\",\n       allow_redirects=false,\n       headers=Dict(\"User-Agent\" => \"Mozilla/5.0 (compatible; MSIE 7.01; Windows NT 5.0)\"))\n\"seq.fasta\"\n\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.get_n_words-Tuple{String,Int64}",
    "page": "Utils",
    "title": "MIToS.Utils.get_n_words",
    "category": "Method",
    "text": "get_n_words{T <: Union{ASCIIString, UTF8String}}(line::T, n::Int) It returns a Vector{T} with the first n (possibles) words/fields (delimited by space, tab or newline). If there is more than n words, the last word returned contains the finals words and the delimiters. The length of the returned vector is n or less (if the number of words is less than n). The last newline character is always removed. This is used for parsing the Stockholm format.\n\njulia> get_n_words(\"#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH\n\", 3)\n3-element Array{String,1}:\n \"#=GR\"\n \"O31698/18-71\"\n \"SS    CCCHHHHHHHHHHHHHHHEEEEEEEEEEEEEEEEHHH\"\n\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.getarray-Tuple{NamedArrays.NamedArray}",
    "page": "Utils",
    "title": "MIToS.Utils.getarray",
    "category": "Method",
    "text": "Getter for the array field of NamedArrays\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.hascoordinates-Tuple{Any}",
    "page": "Utils",
    "title": "MIToS.Utils.hascoordinates",
    "category": "Method",
    "text": "hascoordinates(id) It returns true if id/sequence name has the format: UniProt/start-end (i.e. O83071/192-246)\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.isnotemptyfile-Tuple{Any}",
    "page": "Utils",
    "title": "MIToS.Utils.isnotemptyfile",
    "category": "Method",
    "text": "Returns true if the file exists and isn't empty.\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.lineiterator-Tuple{String}",
    "page": "Utils",
    "title": "MIToS.Utils.lineiterator",
    "category": "Method",
    "text": "Create an iterable object that will yield each line from a stream or string.\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.list2matrix-Tuple{AbstractArray{T,1},Int64}",
    "page": "Utils",
    "title": "MIToS.Utils.list2matrix",
    "category": "Method",
    "text": "Returns a square symmetric matrix from the vector vec. side is the number of rows/columns. The diagonal is not included by default, set to true if there are diagonal elements in the list.\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.matrix2list-Tuple{AbstractArray{T,2}}",
    "page": "Utils",
    "title": "MIToS.Utils.matrix2list",
    "category": "Method",
    "text": "Returns a vector with the part (\"upper\" or \"lower\") of the square matrix mat. The diagonal is not included by default.\n\n\n\n"
},

{
    "location": "Utils_API.html#MIToS.Utils.select_element",
    "page": "Utils",
    "title": "MIToS.Utils.select_element",
    "category": "Function",
    "text": "Selects the first element of the vector. This is useful for unpacking one element vectors. Throws a warning if there are more elements. element_name is element by default, but the name can be changed using the second argument.\n\n\n\n"
},

{
    "location": "Utils_API.html#Methods-and-functions-1",
    "page": "Utils",
    "title": "Methods and functions",
    "category": "section",
    "text": "Modules = [MIToS.Utils]\nPrivate = false\nOrder   = [:function]"
},

]}
