import Base: read

"`FileFormat` is used for write special `parse` (and `read`) methods on it."
abstract type FileFormat end

"""
`download_file` uses **HTTP.jl** instead of system calls to download files
from the web. It takes the file url as first argument and, optionally, a path to save it.
Keyword arguments (ie. `redirect`, `retry`, `readtimeout`)
are are directly passed to to `HTTP.open` (`HTTP.request`).
Use the `headers` keyword argument to pass a `Dict{String,String}` with the
header information.

```jldoctest
julia> using MIToS.Utils

julia> download_file("http://www.uniprot.org/uniprot/P69905.fasta","seq.fasta",
       headers = Dict("User-Agent" =>
                      "Mozilla/5.0 (compatible; MSIE 7.01; Windows NT 5.0)"),
       redirect=true)
"seq.fasta"

```
"""
function download_file(url::AbstractString, filename::AbstractString;
                       headers::Dict{String,String}=Dict{String,String}(),
                       kargs...)
    kargs = _modify_kargs_for_proxy(url; kargs...)
    HTTP.open("GET", url, headers; kargs...) do stream
        open(filename, "w") do fh
           write(fh, stream)
        end
    end
    filename
end

function download_file(url::AbstractString;
                       headers::Dict{String,String}=Dict{String,String}(),
                       kargs...)
    download_file(url, tempname(); headers=headers, kargs...)
end

"""
Helper function that modifies keyword argument to include a proxy,
    the proxy url is taken from the https_proxy and http_proxy enviromental
    variables.
"""
function _modify_kargs_for_proxy(url; kargs...)
    if startswith(lowercase(url), "http://")
        proxy_env_var = "http_proxy"
    elseif startswith(lowercase(url),"https://")
        proxy_env_var = "https_proxy"
    else
        return kargs
    end
    if !(:proxy in keys(kargs)) && proxy_env_var in keys(ENV)
        kw = Dict()
        for (k,v) in kargs
            kw[k] = v
        end
        kw[:proxy] = ENV[proxy_env_var]
        kargs = pairs(kw)
    end
    kargs
end

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
        @warn string(filename, " is empty!")
    end
    filename
end

"Returns `true` if the file exists and isn't empty."
isnotemptyfile(filename) = isfile(filename) && filesize(filename) > 0

# for using with download, since filename doesn't have file extension
function _read(completename::AbstractString,
               filename::AbstractString,
               format::Type{T},
               args...; kargs...) where T <: FileFormat
    check_file(filename)
    if endswith(completename, ".xml.gz") || endswith(completename, ".xml")
        document = parse_file(filename)
        try
            parse(document, T, args...; kargs...)
        finally
            free(document)
        end
    else
        fh = open(filename, "r")
        try
            fh = endswith(completename, ".gz") ? GzipDecompressorStream(fh) : fh
            parse(fh, T, args...; kargs...)
        finally
            close(fh)
        end
    end
end

"""
`read(pathname, FileFormat [, Type [, … ] ] ) -> Type`

This function opens a file in the `pathname` and calls `parse(io, ...)` for
the given `FileFormat` and `Type` on it. If the  `pathname` is an HTTP or FTP URL,
the file is downloaded with `download` in a temporal file.
Gzipped files should end on `.gz`.
"""
function read(completename::AbstractString,
              format::Type{T},
              args...; kargs...) where T <: FileFormat
    if  startswith(completename, "http://")  ||
        startswith(completename, "https://") ||
        startswith(completename, "ftp://")

        filename = download_file(completename)
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
