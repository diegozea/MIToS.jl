"""
`write_file{T<:FileFormat}(filename::AbstractString, object, format::Type{T}, mode::ASCIIString="w")`

This function opens a file with `filename` and `mode` (default: "w")
and writes (`print_file`) the `object` with the given `format`.
Gzipped files should end on `.gz`.
"""
function write_file(
    filename::AbstractString,
    object,
    format::Type{T},
    mode::String = "w",
) where {T<:FileFormat}
    fh = open(filename, mode)
    if endswith(filename, ".gz")
        fh = GzipCompressorStream(fh)
    end
    try
        print_file(fh, object, format)
    finally
        close(fh)
    end
    nothing
end

function Base.write(
    filename::AbstractString,
    object,
    format::Type{T},
    mode::String = "w",
) where {T<:FileFormat}
    Base.depwarn(
        "Using write with $format is deprecated, use write_file instead.",
        :write,
        force = true,
    )
    write_file(filename, object, format, mode)
end

# print_file
# ----------

# Other modules can add their own definition of print_file for their own `FileFormat`s 
# Utils.print_file(io::IO,
print_file(object, format::Type{T}) where {T<:FileFormat} = print_file(stdout, object, T)

function Base.print(fh::IO, object, format::Type{T}) where {T<:FileFormat}
    Base.depwarn(
        "Using print with $format is deprecated, use print_file instead.",
        :print,
        force = true,
    )
    print_file(fh, object, format)
end
