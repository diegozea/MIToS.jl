# Download Pfam
# =============

"""
Download a gzipped stockholm full alignment for the `pfamcode`.
The extension of the downloaded file is `.stockholm.gz` by default.
The `filename` can be changed, but the `.gz` at the end is mandatory.
"""
function downloadpfam(pfamcode::String; filename::String="$pfamcode.stockholm.gz", kargs...)
    @assert endswith(filename,".gz") "filename must end with the .gz extension"
    if occursin(r"^PF\d{5}$"i, pfamcode)
        number = pfamcode[3:end]
        download_file("http://pfam.xfam.org/family/PF$(number)/alignment/full/gzipped",
                      filename; kargs...)
    else
        throw(ErrorException("$pfamcode is not a correct Pfam code"))
    end
end
