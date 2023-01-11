# Download Pfam
# =============

"""
It downloads a gzipped Stockholm alignment from InterPro for the Pfam family 
with the given `pfamcode`. By default, it downloads the `full` Pfam 
alignment. You can use the `alignment` keyword argument to download the 
`seed` or the `uniprot` alignment instead. The extension of the downloaded 
file is `.stockholm.gz` by default; you can change it using the `filename` 
keyword argument, but the `.gz` at the end is mandatory.
"""
function downloadpfam(pfamcode::String; filename::String="$pfamcode.stockholm.gz", 
        alignment::String="full", kargs...)
    if alignment != "full" && alignment != "seed" && alignment != "uniprot"
        throw(ErrorException("alignment must be \"full\", \"seed\" or \"uniprot\""))
    end
    endswith(filename, ".gz") || error("filename must end in .gz")
    if occursin(r"^PF\d{5}$"i, pfamcode)
        download_file("https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/$pfamcode/?annotation=alignment:$alignment&download",
                      filename; kargs...)
    else
        throw(ErrorException("$pfamcode is not a correct Pfam code"))
    end
end
