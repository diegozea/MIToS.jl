# Download Pfam
# =============

"""
Download a gzipped stockholm full alignment for the `pfamcode`.
The extension of the downloaded file is `.stockholm.gz` by default. The `filename` can be changed, but the `.gz` at the end is mandatory.
"""
function downloadpfam(pfamcode::ASCIIString; filename::ASCIIString="$pfamcode.stockholm.gz")
  if ismatch(r"^PF\d{5}$"i, pfamcode) # length(pfamcode)== 7 && ( pfamcode[1:2] == "PF" || pfamcode[1:2] == "pf" )
    #namegz = string(filename, ".gz")
    download(string("http://pfam.xfam.org/family/PF", pfamcode[3:end], "/alignment/full/gzipped"), filename) # namegz)
    #run(`gzip -d $namegz`)
  else
    throw( ErrorException( string(pfamcode, " is not a correct Pfam code") ) )
  end
end


