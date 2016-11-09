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

# AbstractArray Interface
# -----------------------

Base.size(table::ContingencyTable) = size(table.table)

Base.getindex(table::ContingencyTable, i...) = getindex(table.table, i...)

function Base.getindex(table::ContingencyTable, i::Residue)
    index = get_alphabet(table)[i]
    @assert index != 22 "Residue $i isn't in the alphabet"
    getindex(get_table(table), index)
end

function Base.getindex(table::ContingencyTable, i::Residue, j::Residue)
    index_i = get_alphabet(table)[i]
    index_j = get_alphabet(table)[j]
    @assert index_i != 22 "Residue $i isn't in the alphabet"
    @assert index_j != 22 "Residue $j isn't in the alphabet"
    getindex(get_table(table), index_i, index_j)
end

function Base.getindex(table::ContingencyTable, I::Residue...)
    N = length(I)
    indexes = Array(Int,N)
    for i in 1:N
        index = get_alphabet(table)[I[i]]
        @assert index != 22 "Residue $i isn't in the alphabet"
        indexes[i] = index
    end
    getindex(get_table(table), indexes...)
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

function _unsafe_count!{T,N,A}(table::ContingencyTable{T,N,A}, weight::T, i::Int...)
    table.temporal[i...] += weight
end

# Update!
# -------

# Update the table, marginal and total (initialized in 0) using temporal

for (αβ, n) in [(UngappedAlphabet,20), (GappedAlphabet,21)]
    @eval begin

        function _update_table!{T}(table::ContingencyTable{T,1,$αβ})
            array(get_table(table))[:] = table.temporal[1:$n]
            table
        end

        function _update_table!{T}(table::ContingencyTable{T,2,$αβ})
            array(get_table(table))[:] = table.temporal[1:$n, 1:$n]
            table
        end

        function _update_table!{T,N}(table::ContingencyTable{T,N,$αβ})
            array(get_table(table))[:] = table.temporal[(1:$n for i in 1:N)...]
            table
        end

    end
end

_get_indexes(I...) = ((i for i in I)...)

@generated function _update_table!{T,N,A<:ReducedAlphabet}(table::ContingencyTable{T,N,A})
    quote
        temporal = table.temporal
        alphabet = get_alphabet(table)
        freqtable = get_table(table)
        @nloops $N i temporal begin
            indexes = @ncall $N _get_indexes i
            if any(x -> x >= 22, indexes)
                continue
            end
            alphabet_indexes = ((alphabet[x] for x in indexes))
            if any(x -> x >= 22, alphabet_indexes)
                continue
            end
            freqtable[alphabet_indexes...] += @nref $N temporal i
        end
        table
    end
end

function _update_marginals!{T,A}(table::ContingencyTable{T,1,A})
    array(get_marginals(table))[:] = get_table(table)
    table
end

function _update_marginals!{T,A}(table::ContingencyTable{T,2,A})
    sum!(view(array(get_marginals(table)),1,:), array(get_table(table)))
    sum!(reshape(view(array(get_marginals(table)),2,:)
        ,1,length(get_alphabet(table))), array(get_table(table)))
    table
end

function _update_marginals!{T,N,A}(table::ContingencyTable{T,N,A})
    @inbounds for i in 1:N
        array(get_marginals(table))[i,:] =
            sum(get_table(table), ((j for j in 1:N if j != i)...))
    end
    table
end

function _update_total!(table::ContingencyTable)
    @inbounds table.total = sum(view(get_marginals(table), 1, :))
    table
end

function _update!{T,N,A}(table::ContingencyTable{T,N,A})
    _update_table!(table)
    _update_marginals!(table)
    _update_total!(table)
end
