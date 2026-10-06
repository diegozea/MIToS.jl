"""
    print_file(io::IO, msas, format::Type; kwargs...)

Write a collection or iterator of alignments in `Stockholm` or `Clustal` format.
The `msas` iterator must declare an element type other than `Any`, with a `print_file` method
for a single alignment in the requested format. User-defined alignment types use the same
interface as MIToS alignment types.
Alignments are written one at a time, passing any format keywords to each alignment's
writer. An empty typed collection writes nothing. The caller owns the output stream and
any input stream.
"""
function Utils.print_file(io::IO, msas, format::Type{<:Union{Stockholm,Clustal}}; kwargs...)
    Base.IteratorEltype(typeof(msas)) isa Base.HasEltype ||
        throw(ArgumentError("Iterator must declare an alignment `eltype`."))
    eltype(msas) === Any &&
        throw(ArgumentError("Expected an alignment `eltype`; got `Any`."))
    for msa in msas
        # Some scalar values iterate over themselves; stop instead of recursing.
        msa === msas && throw(MethodError(print_file, (io, msas, format)))
        print_file(io, msa, format; kwargs...)
    end
    nothing
end
