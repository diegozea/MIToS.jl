type ContingencyTable{T,N,A} <: AbstractArray{T,N}
    alphabet::A
    temporal::Array{T,N}
    table::NamedArray{T,N,Array{T,N},NTuple{N,OrderedDict{String,Int}}}
    marginals::NamedArray{T,2,Array{T,2},NTuple{2,OrderedDict{String,Int}}}
    total::T
end

# Getters
# -------

@inline getalphabet(table::ContingencyTable) = table.alphabet
@inline gettable(table::ContingencyTable) = table.table
@inline getmarginals(table::ContingencyTable) = table.marginals
@inline gettotal(table::ContingencyTable) = table.total

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
        alphabet = getalphabet(table)
        matrix = array(gettable(table))
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
        alphabet = getalphabet(table)
        matrix = array(gettable(table))
        @nextract $N index d->alphabet[I[d]]
        @_test_index $N index error("There is a Residue outside the alphabet")
        @nref($N, matrix, index) = value
        update_marginals!(table)
    end
end

# Similar
# -------

Base.similar{T,N,A}(table::ContingencyTable{T,N,A}) = ContingencyTable(T, N, table.alphabet)

function Base.similar{T,S,N,A}(table::ContingencyTable{T,N,A}, ::Type{S})
    ContingencyTable(S, N, table.alphabet)
end

# Show
# ====

Base.show(io::IO, ::MIME"text/plain", table::ContingencyTable) = show(io, table)

function Base.show(io::IO, table::ContingencyTable)
    println(io, typeof(table), " : ")
    print(io, "\ntable : ")
    show(io, gettable(table))
    if length(size(table)) != 1
        print(io, "\n\nmarginals : ")
        show(io, getmarginals(table))
    end
    print(io, "\n\ntotal : $(gettotal(table))")
end

# Creation
# --------

@generated function (::Type{ContingencyTable{T,N,A}}){T,N,A}(alphabet::A)
    @assert N > 0 "The dimension should be a natural number"
    quote
        n = length(alphabet)
        # residue_names = names(alphabet)
        # namedict = OrderedDict{String,Int}(residue_names[i] => i for i in 1:n)
        namedict = getnamedict(alphabet)
        dimtable = @ntuple $N k -> n # (n, n, ...)
        dimtemporal = @ntuple $N k -> 22 # (22, 22, ...)
        table = NamedArray(zeros($T, dimtable),
                           @ntuple($N, k -> namedict), # (namedict, namedict, ...)
                           @ntuple($N, k -> "Dim_k")) # ("Dim_1", "Dim_2", ...)
        marginals = NamedArray(zeros($T, $N, n),
                               # OrderedDict{String,Int}("Dim_$i" => i for i in 1:N)
                               (OrderedDict{String,Int}(@ntuple $N k -> "Dim_k"=>k), namedict),
                               ("Dim","Residue"))
        ContingencyTable{T,N,A}(alphabet,
                                zeros(T, dimtemporal),
                                table,
                                marginals,
                                zero(T))
        end
end

function (::Type{ContingencyTable}){T,A}(::Type{T}, N::Int, alphabet::A)
    ContingencyTable{T,N,A}(alphabet)
end

function (::Type{ContingencyTable}){T,N,A}(matrix::AbstractArray{T,N}, alphabet::A)
    n = length(alphabet)
    @assert size(matrix) == ((n for i in 1:N)...) "Matrix size doesn't match alphabet length"
    table = ContingencyTable(T, N, alphabet)
    array(table.table)[:] = matrix
    _update_marginals!(table)
    _update_total!(table)
end

# Update!
# -------

# Update the table, marginal and total using temporal

@generated function _update_table!{T,N,A}(table::ContingencyTable{T,N,A})
    if A <: ReducedAlphabet
        quote
            temporal = table.temporal::Matrix{T}
            alphabet = getalphabet(table)::A
            freqtable = NamedArrays.array(table.table)::Array{T,N}
            @inbounds @nloops $N i temporal begin
                @_test_index $N i continue
                @nextract $N a d->alphabet[i_d]
                @_test_index $N a continue
                @nref($N, freqtable, a) += @nref($N, temporal, i)
            end
            table
        end
    elseif (A === UngappedAlphabet) || (A === GappedAlphabet)
        quote
            temporal = table.temporal::Matrix{T}
            freqtable = NamedArrays.array(table.table)::Array{T,N}
            @inbounds @nloops $N i freqtable begin
                @nref($N, freqtable, i) += @nref($N, temporal, i)
            end
            table
        end
    end
end

@generated function _update_marginals!{T,N,A}(table::ContingencyTable{T,N,A})
    quote
        freqtable = NamedArrays.array(table.table)::Array{T,N}
        marginal  = NamedArrays.array(table.marginals)::Matrix{T}
        @inbounds @nloops $N i freqtable begin
            value = @nref $N freqtable i
            @_marginal($N, marginal, i, value)
        end
        table
    end
end

function _update_total!{T,N,A}(table::ContingencyTable{T,N,A})
    marginals = NamedArrays.array(table.marginals)::Matrix{T}
    ncol = size(marginals,2)
    @inbounds table.total = sum(marginals[1,i] for i in 1:ncol)
    table
end

function update_marginals!{T,N,A}(table::ContingencyTable{T,N,A})
    fill!(NamedArrays.array(table.marginals)::Matrix{T}, zero(T))
    _update_marginals!(table)
    _update_total!(table)
    table
end

function _update!{T,N,A}(table::ContingencyTable{T,N,A})
    _update_table!(table)
    update_marginals!(table)
end

function _cleanup_table!{T,N,A}(table::ContingencyTable{T,N,A})
    fill!(NamedArrays.array(table.table)::Array{T,N}, zero(T))
end

function _cleanup_temporal!{T,N,A}(table::ContingencyTable{T,N,A})
    fill!(table.temporal, zero(T))
end

function cleanup!{T,N,A}(table::ContingencyTable{T,N,A})
    _cleanup_temporal!(table)
    _cleanup_table!(table)
    fill!(NamedArrays.array(table.marginals)::Matrix{T}, zero(T))
    table.total = zero(T)
    table
end

# Fill
# ====

function Base.fill!{T,N,A}(table::ContingencyTable{T,N,A}, value::T)
    fill!(NamedArrays.array(table.table)::Array{T,N}, value)
    update_marginals!(table)
end

Base.fill!{T,N,A}(table::ContingencyTable{T,N,A}, p::AdditiveSmoothing{T}) = fill!(table, p.λ)

@inline Base.fill!{T,N,A}(table::ContingencyTable{T,N,A}, p::NoPseudocount) = table

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

@inline apply_pseudocount!(table::ContingencyTable, p::NoPseudocount) = table

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
`normalize!(p::ContingencyTable)`
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

# Delete Dimensions
# =================

function _list_without_dimensions(len::Int, output_len::Int, dimensions::Int...)
  ndim = length(dimensions)
  @assert (len-ndim) == output_len "$output_len should be = $(len-ndim)"
  index_list = Array(Int, output_len)
  j = 1
  @inbounds for i in 1:len
    if ! (i in dimensions)
      index_list[j] = i
      j += 1
    end
  end
  index_list
end

"""
`delete_dimensions!(out::ContingencyTable, in::ContingencyTable, dimensions::Int...)`

This function fills a ContingencyTable with the counts/probabilities on `in` after the
deletion of `dimensions`. i.e. This is useful for getting Pxy from Pxyz.
"""
function delete_dimensions!{T,N,S,A}(output::ContingencyTable{T,S,A},
                                     input::ContingencyTable{T,N,A},
                                     dimensions::Int...)
  array(output.marginals)[:] = input.marginals[_list_without_dimensions(N, S, dimensions...),:]
  array(output.table)[:] = sum(array(input.table), dimensions)
  output.total = input.total
  output
end

"""
`delete_dimensions(in::ContingencyTable, dimensions::Int...)`

This function creates a ContingencyTable with the counts/probabilities on `in` after the
deletion of `dimensions`. i.e. This is useful for getting Pxy from Pxyz.
"""
function delete_dimensions{T,N,A}(input::ContingencyTable{T,N,A}, dims::Int...)
    delete_dimensions!(ContingencyTable(T, N-length(dims), input.alphabet), input, dims...)
end
