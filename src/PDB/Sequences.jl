# Functions to extract the sequences from PDB structures.
# =======================================================


"""
    is_aminoacid(residue::PDBResidue)
    is_aminoacid(residue_id::PDBResidueIdentifier)

This function returns `true` if the PDB residue is an amino acid residue. It checks if the
residue's three-letter name exists in the `MIToS.Utils.THREE2ONE` dictionary, and
returns `false` otherwise.
"""
is_aminoacid(residue::PDBResidue) = is_aminoacid(residue.id)
is_aminoacid(residue_id::PDBResidueIdentifier) = _is_aminoacid(residue_id.name)
_is_aminoacid(residue_name::String) = haskey(THREE2ONE, residue_name)

function _add_sequence!(chains, key, buf)
    seq = String(take!(buf))
    if !isempty(seq) # do not add empty sequences, for example, if a chain is not selected
        if haskey(chains, key)
            chains[key] *= seq
        else
            chains[key] = seq
        end
    end
end

"""
    modelled_sequences(residue_list::AbstractArray{PDBResidue,N}; 
        model::Union{String,Type{All}}=All, chain::Union{String,Type{All}}=All, 
        group::Union{String,Regex,Type{All}}=All) where N

This function returns an `OrderedDict` where each key is a named tuple (containing the
model and chain identifiers), and each value is the protein sequence corresponding to
the modelled residues in those chains. Therefore, the obtained sequences do not contain
missing residues. All modelled residues are included by default, but those that don't
satisfy specified criteria based on the `model`, `chain`, or `group` keyword arguments
are excluded. One-letter residue names are obtained from the `MIToS.Utils.THREE2ONE`
dictionary for all residue names that return `true` for `is_aminoacid`.
"""
function modelled_sequences(
    residue_list::AbstractArray{PDBResidue,N};
    model::Union{String,Type{All}} = All,
    chain::Union{String,Type{All}} = All,
    group::Union{String,Type{All}} = All,
) where {N}
    chains = OrderedDict{NamedTuple{(:model, :chain),Tuple{String,String}},String}()
    buf = IOBuffer()
    first_residue = first(residue_list)
    key = (model = first_residue.id.model, chain = first_residue.id.chain)
    for res in residue_list
        if !_is(res.id.model, model) ||
           !_is(res.id.chain, chain) ||
           !_is(res.id.group, group) ||
           !is_aminoacid(res)
            continue
        end
        current_key = (model = res.id.model, chain = res.id.chain)
        if current_key != key
            _add_sequence!(chains, key, buf)
            key = current_key
        end
        write(buf, THREE2ONE[res.id.name])
    end
    _add_sequence!(chains, key, buf)
    chains
end
