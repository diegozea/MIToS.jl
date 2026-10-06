"""
    print_file(io::IO, msas, format::Type; kwargs...)

Write a collection or iterator of alignments in `Stockholm` or `Clustal` format.
The `msas` iterator must declare an element type other than `Any`. Each alignment must
be an `AbstractMatrix{Residue}` with a `print_file` method for the requested format.
This also applies to user-defined alignment types.
Alignments are written one at a time, passing any format keywords to each alignment's
writer. An empty typed collection writes nothing. The caller owns the output stream and
any input stream.
"""
function Utils.print_file(io::IO, msas, format::Type{<:Union{Stockholm,Clustal}}; kwargs...)
    # Single alignments use more specific methods; this fallback handles collections.
    # Check the declared element type before consuming the iterator.
    Base.IteratorEltype(typeof(msas)) isa Base.HasEltype ||
        throw(ArgumentError("Iterator must declare an alignment `eltype`."))
    eltype(msas) === Any &&
        throw(ArgumentError("Expected an alignment `eltype`; got `Any`."))
    for msa in msas
        # Require an alignment so unsupported values cannot re-enter this collection method.
        msa isa AbstractMatrix{Residue} || throw(MethodError(print_file, (io, msa, format)))
        # Use each alignment's writer, including methods defined by users.
        print_file(io, msa, format; kwargs...)
    end
    nothing
end
