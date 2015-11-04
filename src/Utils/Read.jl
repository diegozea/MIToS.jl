import Base: read

using GZip
using LightXML

"`Format` is used for write special `parse` (and `read`) methods on it."
abstract Format

"""Returns the `filename`.
Throws an `ErrorException` if the file doesn't exist, or a warning if the file is empty."""
function check_file(filename)
  if !isfile(filename)
    throw(ErrorException(string(filename, " doesn't exist!")))
  elseif filesize(filename) == 0
    warn(string(filename, " is empty!"))
  end
  filename
end

"Returns `true` if the file exists and isn't empty."
isnotemptyfile(filename) = isfile(filename) && filesize(filename) > 0

function _read(completename, filename, format, args...; kargs...) # for using with download, since filename doesn't have file extension
  check_file(filename)
  if endswith(completename, ".xml.gz") || endswith(completename, ".xml")
    document = parse_file(filename)
    try
      parse(document, format, args...; kargs...)
    finally
      free(document)
    end
  else
    io = endswith(completename, ".gz") ? gzopen(filename, "r") : open(filename, "r")
    try
      parse(io, format, args...; kargs...)
    finally
      close(io)
    end
  end
end

"""
```read(pathname, Format [, Type [, … ] ] ) -> Type```

This function opens a file in the `pathname` and calls `parse(io, ...)`for the given `Format` and `Type` on it.
If the  `pathname` is an HTTP or FTP URL, the file is downloaded with `download` in a temporal file.
Gzipped files should end on `.gz`.
"""
function read{T<:Format}(completename::AbstractString, format::Type{T}, args...; kargs...)
  if startswith(completename, "http://") || startswith(completename, "https://") || startswith(completename, "ftp://")
    filename = download(completename)
    try
      _read(completename, filename, format, args...; kargs...)
    finally
      rm(filename)
    end
  else
    _read(completename, completename, format, args...; kargs...) # completename and filename are the same
  end
end
