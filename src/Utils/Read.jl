import Base: read

using GZip
using LightXML

"""
`Format` is used for write special `reader` (and `read`) methods on it.
"""
abstract Format

"""
```read(pathname, Format [, Type [, … ] ] ) -> Type```

This function opens a file in the `pathname` and calls `parse(io, ...)`for the given `Format` and `Type` on it.
If the  `pathname` is an HTTP or FTP URL, the file is downloaded with `download` in a temporal file.
Gzipped files should end on `.gz`.
"""
function read{T<:Format}(pathname::AbstractString, format::Type{T}, args...; kargs...)
  if startswith(pathname, "http://") || startswith(pathname, "https://") || startswith(pathname, "ftp://")
    filename = download(pathname)
    io = endswith(pathname, ".gz") ? gzopen(filename, "r") : open(filename, "r")
    try
      parse(io, format, args...; kargs...)
    finally
      close(io)
      rm(filename)
    end
  elseif endswith(pathname, ".xml.gz") || endswith(pathname, ".xml")
    root = parse_file(pathname) # LightXML.XMLDocument
    parse(root, format, args...; kargs...)
  else
    io = endswith(pathname, ".gz") ? gzopen(pathname, "r") : open(pathname, "r")
    try
      parse(io, format, args...; kargs...)
    finally
      close(io)
    end
  end
end
