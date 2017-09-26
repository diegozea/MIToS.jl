import Base: read

"`Format` is used for write special `parse` (and `read`) methods on it."
abstract Format

"""
`download_file` uses **Requests.jl** instead of system calls to download files from the web.
It takes the file url as first argument and, optionally, a path to save it.
Keyword arguments (ie. `allow_redirects`, `max_redirects`, `timeout`, `headers`)
are are directly passed to to `Requests.get_streaming`.

```julia
julia> download_file("http://www.uniprot.org/uniprot/P69905.fasta","seq.fasta",
       allow_redirects=false,
       headers=Dict("User-Agent" => "Mozilla/5.0 (compatible; MSIE 7.01; Windows NT 5.0)"))
"seq.fasta"

```
"""
function download_file(url::AbstractString, filename::AbstractString; kargs...)
    response = get(url; kargs...)
    save(response, filename)
    filename
end

download_file(url::AbstractString; kargs...) = download_file(url, tempname(); kargs...)

"Create an iterable object that will yield each line from a stream **or string**."
lineiterator(string::String) = eachline(IOBuffer(string))
lineiterator(stream::IO)     = eachline(stream)

"""
Returns the `filename`.
Throws an `ErrorException` if the file doesn't exist, or a warning if the file is empty.
"""
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

# for using with download, since filename doesn't have file extension
function _read(completename, filename, format, args...; kargs...)
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
`read(pathname, Format [, Type [, … ] ] ) -> Type`

This function opens a file in the `pathname` and calls `parse(io, ...)` for
the given `Format` and `Type` on it. If the  `pathname` is an HTTP or FTP URL,
the file is downloaded with `download` in a temporal file.
Gzipped files should end on `.gz`.
"""
function read{T<:Format}(completename::AbstractString, format::Type{T}, args...; kargs...)
    if  startswith(completename, "http://")  ||
        startswith(completename, "https://") ||
        startswith(completename, "ftp://")

        filename = download_file(completename)
        try
            _read(completename, filename, format, args...; kargs...)
        finally
            rm(filename)
        end
    else
        _read(completename, completename, format, args...; kargs...)
            # completename and filename are the same
    end
end
