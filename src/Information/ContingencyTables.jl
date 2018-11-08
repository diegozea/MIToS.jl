# Contingency Tables
# ==================

"""
A `ContingencyTable` is a multidimensional array. It stores the contingency matrix, its
marginal values and total. The type also has an internal and private temporal array and an
alphabet object. It's a parametric type, taking three ordered parameters:

- `T` : The element type of the multidimensional array.
- `N` : It's the dimension of the array and should be an `Int`.
- `A` : This should be a type, subtype of `ResidueAlphabet`, i.e.: `UngappedAlphabet`,
`GappedAlphabet` or `ReducedAlphabet`.

A `ContingencyTable` can be created from an alphabet if all the parameters are given.
Otherwise, you need to give a type, a number (`Val`) and an alphabet. You can also create a
`ContingencyTable` using a matrix and a alphabet. For example:

```julia
ContingencyTable{Float64, 2, UngappedAlphabet}(UngappedAlphabet())
ContingencyTable(Float64, Val{2}, UngappedAlphabet())
ContingencyTable(zeros(Float64,20,20), UngappedAlphabet())
```
"""
mutable struct ContingencyTable{T,N,A} <: AbstractArray{T,N}
    alphabet::A
    temporal::Array{T,N}
    table::NamedArray{T,N,Array{T,N},NTuple{N,OrderedDict{String,Int}}}
    marginals::NamedArray{T,2,Array{T,2},NTuple{2,OrderedDict{String,Int}}}
    total::T
end

# Probability and Counts
# ----------------------

"""
A `Probabilities` object wraps a `ContingencyTable` storing probabilities. It doesn't
perform any check. If the total isn't one, you must use `normalize` or `normalize!`on the
`ContingencyTable` before wrapping it to make the sum of the probabilities equal to one.
"""
mutable struct Probabilities{T,N,A} <: AbstractArray{T,N}
    table::ContingencyTable{T,N,A}
end

"""
A `Counts` object wraps a `ContingencyTable` storing counts/frequencies.
"""
mutable struct Counts{T,N,A} <: AbstractArray{T,N}
    table::ContingencyTable{T,N,A}
end

# Getters

"""
`getcontingencytable` allows to access the wrapped `ContingencyTable` in a `Probabilities`
or `Counts` object.
"""
@inline getcontingencytable(p::Probabilities{T,N,A}) where {T,N,A} = p.table
@inline getcontingencytable(n::Counts{T,N,A}) where {T,N,A} = n.table

for f in (:getalphabet, :gettable, :getmarginals, :gettotal,
          :gettablearray, :getmarginalsarray)
    @eval $(f)(p::Probabilities{T,N,A}) where {T,N,A} = $(f)(getcontingencytable(p))
    @eval $(f)(n::Counts{T,N,A}) where {T,N,A} = $(f)(getcontingencytable(n))
end

# AbstractArray

for f in (:size, :getindex, :setindex!)
    @eval Base.$(f)(p::Probabilities{T,N,A},args...) where {T,N,A} = $(f)(getcontingencytable(p),args...)
    @eval Base.$(f)(n::Counts{T,N,A},args...) where {T,N,A} = $(f)(getcontingencytable(n),args...)
end

# Getters
# -------

"`getalphabet` allows to access the stored alphabet object."
@inline getalphabet(table::ContingencyTable) = table.alphabet

"`gettable` allows to access the table (`NamedArray`)."
@inline gettable(table::ContingencyTable) = table.table

"`getmarginals` allows to access the array with the marginal values (`NamedArray`)."
@inline getmarginals(table::ContingencyTable) = table.marginals

"`gettotal` allows to access the stored total value."
@inline gettotal(table::ContingencyTable) = table.total

Base.sum(table::ContingencyTable) = gettotal(table)

"`gettablearray` allows to access the table (`Array` without names)."
@inline gettablearray(table::ContingencyTable) = getarray(table.table)

"""
`getmarginalsarray` allows to access the array with the marginal values
(`Array` without names).
"""
@inline getmarginalsarray(table::ContingencyTable) = getarray(table.marginals)

# Cartesian (helper functions)
# ----------------------------

"""
`_marginal(1,:A,:i,:value)` generates the expression: `A[i_1, 1] += value`
"""
function _marginal(N::Int, marginal::Symbol, index::Symbol, value::Symbol)
    aexprs = [Expr(:escape, Expr(:(+=), Expr(:ref, marginal, Symbol(index,'_',i), i), :($value))) for i = 1:N]
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
        matrix = getarray(gettable(table))
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
        matrix = getarray(gettable(table))
        @nextract $N index d->alphabet[I[d]]
        @_test_index $N index error("There is a Residue outside the alphabet")
        @nref($N, matrix, index) = value
        update_marginals!(table)
    end
end

# Similar
# -------

Base.similar(table::ContingencyTable{T,N,A}) where {T,N,A} = ContingencyTable(T, Val{N}, table.alphabet)

function Base.similar(table::ContingencyTable{T,N,A}, ::Type{S}) where {T,S,N,A}
    ContingencyTable(S, Val{N}, table.alphabet)
end

# Show
# ====

function Base.show(io::IO, ::MIME"text/plain",
                   table::Union{ContingencyTable{T,N,A},Probabilities{T,N,A},Counts{T,N,A}}) where {T,N,A}
    show(io, table)
end

function Base.show(io::IO, table::ContingencyTable{T,N,A}) where {T,N,A}
    println(io, typeof(table), " : ")
    print(io, "\ntable : ")
    show(io, gettable(table))
    if length(size(table)) != 1
        print(io, "\n\nmarginals : ")
        show(io, getmarginals(table))
    end
    print(io, "\n\ntotal : $(gettotal(table))")
end

function Base.show(io::IO, table::Union{Probabilities{T,N,A},Counts{T,N,A}}) where {T,N,A}
    print(io, typeof(table), " wrapping a ")
    show(io, getcontingencytable(table))
end

# Creation
# --------

@generated function (::Type{ContingencyTable{T,N,A}})(alphabet::A) where {T,N,A}
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
        marginals = NamedArray(zeros($T, n, $N),
                               # OrderedDict{String,Int}("Dim_$i" => i for i in 1:N)
                               (namedict, OrderedDict{String,Int}(@ntuple $N k -> "Dim_k"=>k)),
                               ("Residue","Dim"))
        ContingencyTable{T,N,A}(alphabet,
                                zeros(T, dimtemporal),
                                table,
                                marginals,
                                zero(T))
        end
end

function ContingencyTable(::Type{T}, ::Type{Val{N}}, alphabet::A) where {T,A,N}
    ContingencyTable{T,N,A}(alphabet)
end

function ContingencyTable(matrix::AbstractArray{T,N}, alphabet::A) where {T,N,A}
    n = length(alphabet)
    @assert size(matrix) == ((n for i in 1:N)...,) "Matrix size doesn't match alphabet length"
    table = ContingencyTable(T, Val{N}, alphabet)
    getarray(table.table)[:] = matrix
    _update_marginals!(table)
    _update_total!(table)
end

# Update!
# -------

# Update the table, marginal and total using temporal

@generated function _update_table!(table::ContingencyTable{T,N,A}) where {T,N,A}
    if A <: ReducedAlphabet
        quote
            temporal = table.temporal::Array{T,N}
            alphabet = getalphabet(table)::A
            freqtable = getarray(table.table)::Array{T,N}
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
            temporal = table.temporal::Array{T,N}
            freqtable = getarray(table.table)::Array{T,N}
            @inbounds @nloops $N i freqtable begin
                @nref($N, freqtable, i) += @nref($N, temporal, i)
            end
            table
        end
    end
end

@generated function _update_marginals!(table::ContingencyTable{T,N,A}) where {T,N,A}
        quote
            freqtable = getarray(table.table)::Array{T,N}
            marginal  = getarray(table.marginals)::Matrix{T}
            @inbounds @nloops $N i freqtable begin
                value = @nref $N freqtable i
                @_marginal($N, marginal, i, value)
            end
            table
        end
end

function _update_total!(table::ContingencyTable{T,N,A}) where {T,N,A}
    marginals = getarray(table.marginals)::Matrix{T}
    n = size(marginals,1)
    total = zero(T)
    @inbounds @simd for i in 1:n
        total += marginals[i] # marginals[i,1]
    end
    table.total = total
    table
end

"`update_marginals!` updates the marginal and total values using the table."
function update_marginals!(table::ContingencyTable{T,N,A}) where {T,N,A}
    fill!(getarray(table.marginals)::Matrix{T}, zero(T))
    _update_marginals!(table)
    _update_total!(table)
    table
end

function _update!(table::ContingencyTable{T,N,A}) where {T,N,A}
    _update_table!(table)
    update_marginals!(table)
end

function _cleanup_table!(table::ContingencyTable{T,N,A}) where {T,N,A}
    fill!(getarray(table.table)::Array{T,N}, zero(T))
end

function _cleanup_temporal!(table::ContingencyTable{T,N,A}) where {T,N,A}
    fill!(table.temporal, zero(T))
end

"""
`cleanup!` fills the temporal, table and marginals arrays with zeros.
It also sets total to zero.
"""
function cleanup!(table::ContingencyTable{T,N,A}) where {T,N,A}
    _cleanup_temporal!(table)
    _cleanup_table!(table)
    fill!(getarray(table.marginals)::Matrix{T}, zero(T))
    table.total = zero(T)
    table
end

# Fill
# ====

function Base.fill!(table::ContingencyTable{T,N,A}, value::T) where {T,N,A}
    fill!(getarray(table.table)::Array{T,N}, value)
    update_marginals!(table)
end

Base.fill!(table::ContingencyTable{T,N,A}, p::AdditiveSmoothing{T}) where {T,N,A} = fill!(table, p.λ)

@inline Base.fill!(table::ContingencyTable{T,N,A}, p::NoPseudocount) where {T,N,A} = table

# Apply pseudocount
# =================

# This is faster than array[:] += value
function _sum!(matrix::NamedArray, value)
    matrix_array = getarray(matrix)
    @inbounds for i in eachindex(matrix_array)
        matrix_array[i] += value
    end
    matrix
end

"It adds the `pseudocount` value to the table cells."
function apply_pseudocount!(table::ContingencyTable{T,N,A}, pseudocount::T) where {T,N,A}
    _sum!(table.table, pseudocount)
    update_marginals!(table)
end

function apply_pseudocount!(table::ContingencyTable{T,N,A}, p::AdditiveSmoothing{T}) where {T,N,A}
    apply_pseudocount!(table, p.λ)
end

@inline apply_pseudocount!(table::ContingencyTable, p::NoPseudocount) = table

# Normalize
# =========

# This is faster than array[:] /= value
function _div!(matrix::NamedArray, value)
    matrix_array = getarray(matrix)
    @inbounds for i in eachindex(matrix_array)
        matrix_array[i] /= value
    end
    matrix
end

"`normalize!` makes the sum of the frequencies to be one, in place."
function LinearAlgebra.normalize!(table::ContingencyTable{T,N,A}) where {T,N,A}
    if table.total != one(T)
        _div!(table.table, table.total)
        update_marginals!(table)
    end
    table
end

"`normalize` returns another table where the sum of the frequencies is one."
LinearAlgebra.normalize(table::ContingencyTable{T,N,A}) where {T,N,A} = normalize!(deepcopy(table))

# Delete Dimensions
# =================

function _list_without_dimensions(len::Int, output_len::Int, dimensions::Int...)
  ndim = length(dimensions)
  @assert (len-ndim) == output_len "$output_len should be = $(len-ndim)"
  index_list = Array{Int}(undef, output_len)
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
function delete_dimensions!(output::ContingencyTable{T,S,A},
                            input::ContingencyTable{T,N,A},
                            dimensions::Int...) where {T,N,S,A}
  output_marginals = getarray(output.marginals)
  output_table = getarray(output.table)
  input_marginals = getarray(input.marginals)
  input_table = getarray(input.table)
  output_marginals[:] = input_marginals[:,_list_without_dimensions(N, S, dimensions...)]
  output_table[:] = sum(input_table, dims=dimensions)
  output.total = input.total
  output
end

"""
`delete_dimensions(in::ContingencyTable, dimensions::Int...)`

This function creates a ContingencyTable with the counts/probabilities on `in` after the
deletion of `dimensions`. i.e. This is useful for getting Pxy from Pxyz.
"""
function delete_dimensions(input::ContingencyTable{T,N,A}, dims::Vararg{Int,I}) where {T,N,A,I}
    delete_dimensions!(ContingencyTable(T, Val{N-I}, input.alphabet), input, dims...)
end

for tp in (:Probabilities, :Counts)
    @eval begin
        function delete_dimensions!(output::$(tp){T,S,A}, input::$(tp){T,N,A},
                                    dimensions::Int...) where {T,N,S,A}
            delete_dimensions!(getcontingencytable(output), getcontingencytable(input), dimensions...)
            output
        end

        function delete_dimensions(input::$(tp){T,N,A}, dims::Vararg{Int,I}) where {T,N,A,I}
            output = delete_dimensions(getcontingencytable(input), dims...)
            $(tp)(output)
        end
    end
end
