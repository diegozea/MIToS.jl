import Base: read

using GZip

"""
`Format` is used for write special `reader` (and `read`) methods on it.
"""
abstract Format

reader(io::IO) = throw(ErrorException("You need to indicate the format of the file (and maybe the output type)"))

"""
```read(pathname, Format [, Type [, … ] ] ) -> Type```

This function opens a file in the `pathname` and calls `reader`for the given `Format` and `Type` on it.
If the  `pathname` is an HTTP or FTP URL, the file is downloaded with `download` in a temporal file.
Gzipped files should end on `.gz`.
"""
function read(pathname::AbstractString, args...; kargs...)
  if startswith(pathname, "http://") || startswith(pathname, "https://") || startswith(pathname, "ftp://")
    filename = download(pathname)
    io = endswith(pathname, ".gz") ? gzopen(filename, "r") : open(filename, "r")
    try
      reader(io, args...; kargs...)
    finally
      close(io)
      rm(filename)
    end
  else
    io = endswith(pathname, ".gz") ? gzopen(pathname, "r") : open(pathname, "r")
    try
      reader(io, args...; kargs...)
    finally
      close(io)
    end
  end
end
