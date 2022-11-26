# Download Pfam
# =============

"""
Download a gzipped stockholm full alignment for the `pfamcode`.
The extension of the downloaded file is `.stockholm.gz` by default.
The `filename` can be changed, but the `.gz` at the end is mandatory.
"""
function downloadpfam(pfamcode::String; filename::String="$pfamcode.alignment.uniprot.gz", kargs...)
    m = match(r"alignment\.([a-z]*)\.gz$", filename)
    m === nothing && error("filename must end in alignment.(seed|full|uniprot).gz")
    mode = m.captures[1]
    if occursin(r"^PF\d{5}$"i, pfamcode)
        download_file("https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/$pfamcode/?annotation=alignment:$mode&download",
                      filename; kargs...)
    else
        throw(ErrorException("$pfamcode is not a correct Pfam code"))
    end
end
