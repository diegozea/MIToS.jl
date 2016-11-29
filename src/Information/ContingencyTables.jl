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

# Update the table, marginal and total (initialized in 0) using temporal

for (αβ, n) in [(UngappedAlphabet,20), (GappedAlphabet,21)]
    @eval begin

        function _update_table!{T}(table::ContingencyTable{T,1,$αβ})
            array(table.table)[:] = table.temporal[1:$n]
            table
        end

        function _update_table!{T}(table::ContingencyTable{T,2,$αβ})
            array(table.table)[:] = table.temporal[1:$n, 1:$n]
            table
        end

        function _update_table!{T,N}(table::ContingencyTable{T,N,$αβ})
            array(table.table)[:] = table.temporal[(1:$n for i in 1:N)...]
            table
        end

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

function _update!{T,N,A}(table::ContingencyTable{T,N,A})
    _update_table!(table)
    _update_marginals!(table)
    _update_total!(table)
    table
end

function update_marginals!{T,N,A}(table::ContingencyTable{T,N,A})
    fill!(table.marginals, zero(T))
    _update_marginals!(table)
    _update_total!(table)
    table
end

# Fill
# ====

function Base.fill!{T,N,A}(table::ContingencyTable{T,N,A}, value::T)
    fill!(table.table, value)
    update_marginals!(table)
end

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
The marginals are updated in the normalization.
"""
function Base.normalize!{T,N,A}(table::ContingencyTable{T,N,A})
    if table.total != T(1.0)
        _div!(table.table, table.total)
    end
    update_marginals!(table)
end

