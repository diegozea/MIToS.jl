"""
Validate an object before opening its output file. Formats can specialize this check.
"""
_validate_write(object, format::Type{<:FileFormat}) = nothing

"""
    write_file(filename::AbstractString, object, format::Type, mode::String = "w")

This function opens a file with `filename` and `mode` (default: "w")
and writes (`print_file`) the `object` with the given `format`.
Gzipped files should end on `.gz`.

For `Stockholm` and `Clustal`, `object` can also be a collection or iterator with a
declared alignment element type (`eltype(object) <: AbstractMatrix{Residue}`). Alignments
are written one at a time. Unknown or incompatible element types raise an `ArgumentError`
before the file is opened. An empty typed collection writes an empty file.
"""
function write_file(
    filename::AbstractString,
    object,
    format::Type{T},
    mode::String = "w",
) where {T<:FileFormat}
    _validate_write(object, format)
    fh = open(filename, mode)
    try
        if endswith(filename, ".gz")
            fh = GzipCompressorStream(fh)
            write(fh, "") # Start a valid gzip stream even when there are no alignments.
        end
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
