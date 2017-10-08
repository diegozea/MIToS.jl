import Base: write

"""
`write{T<:Format}(filename::AbstractString, object, format::Type{T}, mode::ASCIIString="w")`

This function opens a file with `filename` and `mode` (default: "w")
and writes (`print`) the `object` with the given `format`.
Gzipped files should end on `.gz`.
"""
function write(filename::AbstractString, object, format::Type{T},
               mode::String="w") where T<:Format
    fh = endswith(filename, ".gz") ? GZip.open(filename, mode) : open(filename, mode)
    try
        print(fh, object, format)
    finally
        close(fh)
    end
    nothing
end
