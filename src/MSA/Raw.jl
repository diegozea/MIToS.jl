struct Raw <: Format end

# Raw Parser
# ==========

function _get_seqs(io::Union{IO, AbstractString})
    SEQS = String[]
    for line in lineiterator(io)
        push!(SEQS, chomp(line))
    end
    SEQS
end

function Base.parse{T <: MSAMatrix}(io::Union{IO, AbstractString},
                                    format::Type{Raw},
                                    output::Type{T};
                                    deletefullgaps::Bool=true)::T
    SEQS = _get_seqs(io)
    _strings_to_msa(T, SEQS, deletefullgaps)
end

function Base.parse(io::Union{IO, AbstractString},
                    format::Type{Raw},
                    output::Type{AnnotatedMultipleSequenceAlignment};
                    generatemapping::Bool=false,
                    deletefullgaps::Bool=true,
                    keepinserts::Bool=false)
    SEQS = _get_seqs(io)
    annot = Annotations()
    if keepinserts
        _keepinserts!(SEQS, annot)
    end
    if generatemapping
        MSA, MAP = _to_msa_mapping(SEQS) # It checks sequence lengths
        setannotfile!(annot, "NCol", string(size(MSA,2)))
        setannotfile!(annot, "ColMap", join(vcat(1:size(MSA,2)), ','))
        for i in 1:length(SEQS)
            setannotsequence!(annot, string(i), "SeqMap", MAP[i])
        end
    else
        MSA = NamedArray(convert(Matrix{Residue}, SEQS)) # It checks sequence lengths
    end
    msa = AnnotatedMultipleSequenceAlignment(MSA, annot)
    if deletefullgaps
        deletefullgapcolumns!(msa)
    end
    msa
end

function Base.parse(io::Union{IO, AbstractString},
                    format::Type{Raw},
                    output::Type{MultipleSequenceAlignment};
                    deletefullgaps::Bool=true)
    msa = parse(io, Raw, NamedArray{Residue,2}; deletefullgaps=deletefullgaps)
    MultipleSequenceAlignment(msa)
end

function Base.parse(io::Union{IO, AbstractString}, format::Type{Raw};
                    deletefullgaps::Bool=true)
    parse(io, Raw, AnnotatedMultipleSequenceAlignment; deletefullgaps=deletefullgaps)
end

# Print Raw
# =========

function Base.print(io::IO, msa::AbstractMatrix{Residue}, format::Type{Raw})
    for i in 1:nsequences(msa)
        println(io, stringsequence(msa, i))
    end
end

Base.print(msa::AbstractMatrix{Residue}, format::Type{Raw}) = print(STDOUT, msa, Raw)
