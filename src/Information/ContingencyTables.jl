type ContingencyTable{T,N,A} <: AbstractArray{T,N}
    alphabet::ResidueAlphabet
    temporal::Array{T,N}
    table::NamedArray{T,N}
    marginals::NamedArray{T,2}
    total::T
end

# Getters
# -------

@inline get_alphabet(table::ContingencyTable) = table.alphabet
@inline get_table(table::ContingencyTable) = table.table
@inline get_marginals(table::ContingencyTable) = table.marginals
@inline get_total(table::ContingencyTable) = table.total

# Cartesian (helper functions)
# ----------------------------

"""
`_marginal(1,:A,:i,:value)` generates the expression: `A[1,i_1] += value`
"""
function _marginal(N::Int, marginal::Symbol, index::Symbol, value::Symbol)
    aexprs = [Expr(:escape, Expr(:(+=), Expr(:ref, marginal, i, Symbol(index,'_',i)), :($value))) for i = 1:N]
    Expr(:block, aexprs...)
end

macro _marginal(N, marginal, index, value)
    _marginal(N, marginal, index, value)
end

"""
`_test_index(1, i, continue)` generates the expression: `i_1 >= 22 && continue`
"""
function _test_index(N::Int, index::Symbol, expr::Expr)
    aexprs = [Expr(:escape, :($(Symbol(index,"_",i)) >= 22 && $expr)) for i = 1:N]
    Expr(:block, aexprs...)
end

macro _test_index(N, index, expr)
    _test_index(N, index, expr)
end

# AbstractArray Interface
# -----------------------

Base.size(table::ContingencyTable) = size(table.table)

Base.getindex(table::ContingencyTable, i...) = getindex(table.table, i...)

@generated function Base.getindex(table::ContingencyTable, I::Residue...)
    N = length(I)
    quote
        alphabet = get_alphabet(table)
        matrix = array(get_table(table))
        # index_1 = alphabet[I[1]]
        # index_2 ...
        @nextract $N index d->alphabet[I[d]]
        # index_1 >= 22 error("There is a Residue outside the alphabet")
        # index_2 ...
        @_test_index $N index error("There is a Residue outside the alphabet")
        # getindex(matrix, index_1, index_2...
        @nref $N matrix index
    end
end

function Base.setindex!(table::ContingencyTable, value, i...)
    setindex!(table.table, value, i...)
    update_marginals!(table)
end

@generated function Base.setindex!(table::ContingencyTable, value, I::Residue...)
    N = length(I)
    quote
        alphabet = get_alphabet(table)
        matrix = array(get_table(table))
        @nextract $N index d->alphabet[I[d]]
        @_test_index $N index error("There is a Residue outside the alphabet")
        @nref($N, matrix, index) = value
        update_marginals!(table)
    end
end

# Show
# -----

Base.show(io::IO, ::MIME"text/plain", table::ContingencyTable) = show(io, table)

function Base.show(io::IO, table::ContingencyTable)
    println(io, typeof(table), " : ")
    print(io, "\ntable : ")
    show(io, get_table(table))
    if length(size(table)) != 1
        print(io, "\n\nmarginals : ")
        show(io, get_marginals(table))
    end
    print(io, "\n\ntotal : $(get_total(table))")
end

# Creation
# --------

function (::Type{ContingencyTable}){T,A}(::Type{T}, N::Int, alphabet::A)
    @assert N > 0 "The dimension should be a natural number"
    n = length(alphabet)
    residue_names = names(alphabet)
    dim = ((n for i in 1:N)...)
    table = NamedArray(zeros(T, dim))
    for i in 1:N
        setnames!(table, residue_names, i)
    end
    setdimnames!(table, (("Dim_$i" for i in 1:N)...))
    marginals = NamedArray(zeros(T, N, n))
    setnames!(marginals, residue_names, 2)
    setdimnames!(marginals, ("Dim","Residue"))
    ContingencyTable{T,N,A}(alphabet,
                            zeros(T, (22 for i in 1:N)...),
                            table,
                            marginals,
                            zero(T))
end

# Unsafe count
# ------------

# Count into the temporal field, without updating table, marginals and total.

@inline function _unsafe_count!{T,N,A}(table::ContingencyTable{T,N,A}, weight::T, i::Int...)
    @inbounds table.temporal[i...] += weight
end

# Update!
# -------

# Update the table, marginal and total using temporal

@generated function _update_table!{T,N,A<:Union{UngappedAlphabet,GappedAlphabet}}(table::ContingencyTable{T,N,A})
    quote
        temporal = table.temporal
        freqtable = array(table.table)
        @inbounds @nloops $N i freqtable begin
            @nref($N, freqtable, i) += @nref($N, temporal, i)
        end
        table
    end
end

@generated function _update_table!{T,N,A<:ReducedAlphabet}(table::ContingencyTable{T,N,A})
    quote
        temporal = table.temporal
        alphabet = get_alphabet(table)
        freqtable = array(table.table)
        @inbounds @nloops $N i temporal begin
            @_test_index $N i continue
            @nextract $N a d->alphabet[i_d]
            @_test_index $N a continue
            @nref($N, freqtable, a) += @nref($N, temporal, i)
        end
        table
    end
end

@generated function _update_marginals!{T,N,A}(table::ContingencyTable{T,N,A})
    quote
        freqtable = array(table.table)
        marginal  = array(table.marginals)
        @inbounds @nloops $N i freqtable begin
            value = @nref $N freqtable i
            @_marginal($N, marginal, i, value)
        end
        table
    end
end

function _update_total!(table::ContingencyTable)
    @inbounds table.total = sum(view(table.marginals, 1, :))
    table
end

function update_marginals!{T,N,A}(table::ContingencyTable{T,N,A})
    fill!(table.marginals, zero(T))
    _update_marginals!(table)
    _update_total!(table)
    table
end

function _update!{T,N,A}(table::ContingencyTable{T,N,A})
    _update_table!(table)
    update_marginals!(table)
end

_cleanup_temporal!{T,N,A}(table::ContingencyTable{T,N,A}) = fill!(table.temporal, zero(T))

# Fill
# ====

function Base.fill!{T,N,A}(table::ContingencyTable{T,N,A}, value::T)
    fill!(table.table, value)
    update_marginals!(table)
end

Base.fill!{T,N,A}(table::ContingencyTable{T,N,A}, p::AdditiveSmoothing{T}) = fill!(table, p.λ)

# Apply pseudocount
# =================

# This is faster than array[:] += value
function _sum!(matrix::NamedArray, value)
    matrix_array = array(matrix)
    @inbounds for i in eachindex(matrix_array)
        matrix_array[i] += value
    end
    matrix
end

function apply_pseudocount!{T,N,A}(table::ContingencyTable{T,N,A}, pseudocount::T)
    _sum!(table.table, pseudocount)
    update_marginals!(table)
end

function apply_pseudocount!{T,N,A}(table::ContingencyTable{T,N,A}, p::AdditiveSmoothing{T})
    apply_pseudocount!(table, p.λ)
end

# Normalize
# =========

# This is faster than array[:] /= value
function _div!(matrix::NamedArray, value)
    matrix_array = array(matrix)
    @inbounds for i in eachindex(matrix_array)
        matrix_array[i] /= value
    end
    matrix
end

"""
`normalize!(p::ResidueCount)`
This function makes the sum of the frequencies to be one.
"""
function Base.normalize!{T,N,A}(table::ContingencyTable{T,N,A})
    if table.total != one(T)
        _div!(table.table, table.total)
        _div!(table.marginals, table.total)
        table.total = one(T)
    end
    table
end

# Counting
# ======

@generated function _temporal_counts!{T,N,A}(counts::ContingencyTable{T,N,A}, weights,
                                  seqs::AbstractVector{Residue}...)
    quote
        @assert N == length(seqs) "Number of residue arrays and table dimension doesn't match."
        # seq_1 = seqs[1]
        # seq_2 = ...
        @nextract $N seq d -> seqs[d]
        len = length(seq_1)
        # @assert len == length(seq_1) "Residue arrays have different lengths"
        # @assert len == length(seq_2) ...
        @nexprs $N d -> @assert len == length(seq_d) "Residue arrays have different lengths."
        if isa(weights, AbstractArray)
            @assert len == length(weights) "Residue array and weights sizes doesn't match."
        end
        _cleanup_temporal!(counts)
        temporal = counts.temporal
        @inbounds @simd for index in 1:length(seq_1)
            # temporal[Int(seq_1[index]), Int(seq_2... += getweight(weights, index)
            @nref($N, temporal, d -> Int(seq_d[index])) += getweight(weights, index)
        end
        counts
    end
end

function count!{T,N,A}(table::ContingencyTable{T,N,A}, weights, seqs::AbstractVector{Residue}...)
    _temporal_counts!(table, weights, seqs...)
    _update!(table)
end
