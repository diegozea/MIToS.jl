"""
Check whether a group of multiple sequence alignments can be saved in the requested
file format. Each declared alignment type must have a `print_file` method for that format.
The check does not read the alignments or inspect their sequences. `write_file` runs it
before creating or overwriting the output file.
Single alignments and groups with their own `print_file` method need no additional check.
"""
function Utils._validate_write(
    msas,
    format::Type{F},
    ::Type{S},
) where {F<:Union{Stockholm,Clustal},S<:IO}
    # Only the collection fallback needs this check before write_file opens the file.
    collection_writer = which(print_file, (S, Any, Type{F}))
    which(print_file, (S, typeof(msas), Type{F})) === collection_writer || return nothing
    Base.IteratorEltype(typeof(msas)) isa Base.HasEltype ||
        throw(ArgumentError("Iterator must declare an alignment eltype."))
    T = eltype(msas)
    T === Union{} && return nothing
    all(
        type -> which(print_file, (S, type, Type{F})) !== collection_writer,
        Base.uniontypes(T),
    ) || throw(ArgumentError("Expected an alignment eltype; got $T."))
    nothing
end

"""
    print_file(io::IO, msas, format::Type; kwargs...)

Write a collection or iterator of alignments in `Stockholm` or `Clustal` format.
The iterator must declare an element type other than `Any`, with a `print_file` method
for a single alignment in the requested format. User-defined alignment types use the same
interface as MIToS alignment types.
Alignments are written one at a time, passing any format keywords to each alignment's
writer. An empty typed collection writes nothing. The caller owns the output stream and
any input stream.
"""
function Utils.print_file(io::IO, msas, format::Type{<:Union{Stockholm,Clustal}}; kwargs...)
    Utils._validate_write(msas, format, typeof(io))
    for msa in msas
        print_file(io, msa, format; kwargs...)
    end
    nothing
end
