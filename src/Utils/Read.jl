import Base: read

"""
`FileFormat` is used for defile special `parse_file` (called by `read_file`) and
`print_file` (called by `read_file`) methods for different file formats.
"""
abstract type FileFormat end

"""
Return whether `source` is an HTTP, HTTPS or FTP URL.
"""
_is_url(source::AbstractString) =
    any(prefix -> startswith(source, prefix), ("http://", "https://", "ftp://"))

"""
Create a temporary download filename, preserving a `.gz` suffix.
"""
_download_tempname(url::AbstractString) = tempname() * (endswith(url, ".gz") ? ".gz" : "")

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

julia> download_file(\"https://www.uniprot.org/uniprot/P69905.fasta\", \"seq.fasta\")
"seq.fasta"
```
"""
function download_file(url::AbstractString, filename::AbstractString; kargs...)
    retry(_download_file, delays = ExponentialBackOff(n = 5))(url, filename; kargs...)
end

function download_file(url::AbstractString; kargs...)
    download_file(url, _download_tempname(url); kargs...)
end

"""
Create an iterable object that will yield each line from a stream **or string**.
"""
@inline lineiterator(string::String) = eachline(IOBuffer(string))
@inline lineiterator(stream::IO) = eachline(stream)

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

"""
Wrap `io` in a gzip decompressor when `source` ends in `.gz`; otherwise return `io`.
"""
_input_stream(io::IO, source::AbstractString) =
    endswith(source, ".gz") ? GzipDecompressorStream(io) : io

"""
Close an optional stream and remove its temporary download, even if closing fails.
"""
function _close_input(io::Union{Nothing,IO}, temporary::Union{Nothing,String})
    try
        io === nothing || close(io)
    finally
        temporary === nothing || rm(temporary; force = true)
    end
    nothing
end

function _get_xml_document(filename::AbstractString)
    if endswith(filename, ".gz")
        _check_gzip_file(filename)
        open(filename, "r") do fh
            xml = read(_input_stream(fh, filename), String)
            return LightXML.parse_string(xml)
        end
    else
        return LightXML.parse_file(filename)
    end
end

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
        document = _get_xml_document(filename)
        try
            parse_file(document, T, args...; kargs...)
        finally
            LightXML.free(document)
        end
    else
        open(filename, "r") do fh
            fh = _input_stream(fh, completename)
            _read(fh, T, args...; kargs...)
        end
    end
end

"""
Read an object from an open file. Formats can specialize this method to check for
additional records without changing `parse_file` or how files are opened and closed.
"""
function _read(io::IO, format::Type{<:FileFormat}, args...; kwargs...)
    parse_file(io, format, args...; kwargs...)
end

"""
`read_file(pathname, FileFormat [, Type [, … ] ] ) -> Type`

This function opens a file in the `pathname` and calls `parse_file(io, ...)` for
the given `FileFormat` and `Type` on it. If the  `pathname` is an HTTP or FTP URL,
the file is downloaded with `download` in a temporal file.
Gzipped files should end on `.gz`.

For Stockholm and Clustal files, only the first alignment is returned. A warning is
shown if another alignment is found; use [`eachmsa`](@ref MIToS.MSA.eachmsa) to read them all.
"""
function read_file(
    completename::AbstractString,
    format::Type{T},
    args...;
    kargs...,
) where {T<:FileFormat}
    if _is_url(completename)
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
    Base.depwarn(
        "Using read with $format is deprecated, use read_file instead.",
        :read,
        force = true,
    )
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
    Base.depwarn(
        "Using parse with $format is deprecated, use parse_file instead.",
        :parse,
        force = true,
    )
    parse_file(io, format, args...; kargs...)
end

# A placeholder to define the function name so that other modules can add their own 
# definition of parse_file for their own `FileFormat`s
function parse_file end
