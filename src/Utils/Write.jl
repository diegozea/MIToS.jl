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

"""
    print_file(io::IO, object, format::Type{T}) where T<:FileFormat
    print_file(object, format::Type{T}) where T<:FileFormat

Prints the `object` to the specified `IO` stream (or `stdout` if `io` is omitted)
in the given `FileFormat`.

This is a generic function. Other modules should extend it by defining methods for
their specific `FileFormat` subtypes and the types of objects they handle.

**Example of extending `print_file`:**

```julia
# In another module, after defining MyFormat <: FileFormat and MyObject:
# function MIToS.Utils.print_file(io::IO, data::MyObject, ::Type{MyFormat})
#     # Custom printing logic here, writing to io
#     println(io, "Custom format output for MyObject")
# end
```
"""
function print_file end # Generic function declaration

# Method defaulting to stdout
"""
    print_file(object, format::Type{T}) where T<:FileFormat

Prints the `object` to `stdout` in the given `FileFormat`.
This is a convenience method that calls `print_file(stdout, object, format)`.
"""
print_file(object, format::Type{T}) where {T<:FileFormat} = print_file(stdout, object, T)

# Other modules can add their own definition of print_file for their own `FileFormat`s
# e.g. MIToS.Utils.print_file(io::IO, data::MyType, ::Type{MyFormat})

function Base.print(fh::IO, object, format::Type{T}) where {T<:FileFormat}
    Base.depwarn(
        "Using print with $format is deprecated, use print_file instead.",
        :print,
        force = true,
    )
    print_file(fh, object, format)
end
