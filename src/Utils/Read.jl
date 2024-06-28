import Base: read

"""
`FileFormat` is used for defile special `parse_file` (called by `read_file`) and
`print_file` (called by `read_file`) methods for different file formats.
"""
abstract type FileFormat end

"""
This function raises an error if a GZip file doesn't have the 0x1f8b magic number.
"""
function _check_gzip_file(filename)
    if endswith(filename, ".gz")
        open(filename, "r") do fh
            magic = read(fh, UInt16)
            # 0x1f8b is the magic number for GZip files
            # However, some files use 0x8b1f.
            # For example, the file test/data/18gs.xml.gz uses 0x8b1f.
            if magic != 0x1f8b && magic != 0x8b1f
                throw(ErrorException("$filename is not a GZip file!"))
            end
        end
    end
    filename
end

function _download_file(url::AbstractString, filename::AbstractString; kargs...)
    with_logger(ConsoleLogger(stderr, Logging.Warn)) do
        Downloads.download(url, filename; kargs...)
    end
    _check_gzip_file(filename)
end

"""
`download_file` uses **Downloads.jl** to download files from the web. It takes the file
url as first argument and, optionally, a path to save it.
Keyword arguments are are directly passed to to `Downloads.download`.

```jldoctest
julia> using MIToS.Utils

julia> download_file(
           "http://www.uniprot.org/uniprot/P69905.fasta",
           "seq.fasta",
           headers = Dict("User-Agent" => "Mozilla/5.0 (compatible; MSIE 7.01; Windows NT 5.0)"),
       )
"seq.fasta"
```
"""
function download_file(url::AbstractString, filename::AbstractString; kargs...)
    retry(_download_file, delays = ExponentialBackOff(n = 5))(url, filename; kargs...)
end

function download_file(url::AbstractString; kargs...)
    name = tempname()
    if endswith(url, ".gz")
        name *= ".gz"
    end
    download_file(url, name; kargs...)
end

"""
Create an iterable object that will yield each line from a stream **or string**.
"""
lineiterator(string::String) = eachline(IOBuffer(string))
lineiterator(stream::IO) = eachline(stream)

"""
Returns the `filename`.
Throws an `ErrorException` if the file doesn't exist, or a warning if the file is empty.
"""
function check_file(filename)
    if !isfile(filename)
        throw(ErrorException(string(filename, " doesn't exist!")))
    elseif filesize(filename) == 0
        @warn string(filename, " is empty!")
    end
    filename
end

"""
Returns `true` if the file exists and isn't empty.
"""
isnotemptyfile(filename) = isfile(filename) && filesize(filename) > 0

# for using with download, since filename doesn't have file extension
function _read(
    completename::AbstractString,
    filename::AbstractString,
    format::Type{T},
    args...;
    kargs...,
) where {T<:FileFormat}
    check_file(filename)
    if endswith(completename, ".xml.gz") || endswith(completename, ".xml")
        document = LightXML.parse_file(filename)
        try
            parse_file(document, T, args...; kargs...)
        finally
            LightXML.free(document)
        end
    else
        fh = open(filename, "r")
        try
            fh = endswith(completename, ".gz") ? GzipDecompressorStream(fh) : fh
            parse_file(fh, T, args...; kargs...)
        finally
            close(fh)
        end
    end
end

"""
`read_file(pathname, FileFormat [, Type [, … ] ] ) -> Type`

This function opens a file in the `pathname` and calls `parse_file(io, ...)` for
the given `FileFormat` and `Type` on it. If the  `pathname` is an HTTP or FTP URL,
the file is downloaded with `download` in a temporal file.
Gzipped files should end on `.gz`.
"""
function read_file(
    completename::AbstractString,
    format::Type{T},
    args...;
    kargs...,
) where {T<:FileFormat}
    if startswith(completename, "http://") ||
       startswith(completename, "https://") ||
       startswith(completename, "ftp://")

        filename =
            download_file(completename, headers = Dict("Accept-Encoding" => "identity"))
        try
            _read(completename, filename, T, args...; kargs...)
        finally
            rm(filename)
        end
    else
        # completename and filename are the same
        _read(completename, completename, T, args...; kargs...)
    end
end

function read(
    name::AbstractString,
    format::Type{T},
    args...;
    kargs...,
) where {T<:FileFormat}
    Base.depwarn("Using read with $format is deprecated, use read_file instead.", :read, force=true)
    read_file(name, format, args...; kargs...)
end

# parse_file
# ----------

function Base.parse(
    io::Union{IO,AbstractString},
    format::Type{T},
    args...;
    kargs...,
) where {T<:FileFormat}
    Base.depwarn("Using parse with $format is deprecated, use parse_file instead.", :parse, force=true)
    parse_file(io, format, args...; kargs...)
end

# A placeholder to define the function name so that other modules can add their own 
# definition of parse_file for their own `FileFormat`s
function parse_file end
