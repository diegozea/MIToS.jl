import Base: write

using GZip

"""
```write{T<:Format}(filename::AbstractString, object, format::Type{T}, mode::ASCIIString="w")```

This function opens a file with `filename` and `mode` (default: "w") and writes (`print`) the `object` with the given `format`.
Gzipped files should end on `.gz`.
"""
function write{T<:Format}(filename::AbstractString, object, format::Type{T}, mode::ASCIIString="w")
    fh = endswith(filename, ".gz") ? GZip.open(filename, mode) : open(filename, mode)
    print(fh, object, format)
    close(fh)
    nothing
end
