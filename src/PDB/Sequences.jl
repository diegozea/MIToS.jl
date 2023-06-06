# Functions to extract the sequences from PDB structures.
# =======================================================

function _add_sequence!(buf, chains, key)
    seq = String(take!(buf))
    if haskey(chains, key)
        chains[key] *= seq
    else
        chains[key] = seq
    end
end

"""
    modelled_sequences(residue_list::AbstractArray{PDBResidue,N}; 
        model::Union{String,Type{All}}=All, chain::Union{String,Type{All}}=All, 
        group::Union{String,Regex,Type{All}}=All) where N

This function generates a dictionary where each key-value pair consists of a named tuple 
(representing the `model` and `chain` identifiers) and a corresponding protein sequence. 
All residues are included by default, but those that don't satisfy specified criteria 
based on `model`, `chain`, and `group` keyword arguments are excluded. If a residue name 
isn't located in the `MIToS.Utils.THREE2ONE` dictionary, the function triggers a warning.
"""
function modelled_sequences(residue_list::AbstractArray{PDBResidue,N}; 
        model::Union{String,Type{All}}=All, chain::Union{String,Type{All}}=All, group::Union{String,Type{All}}=All) where N
    chains = OrderedDict{NamedTuple{(:model, :chain), Tuple{String, String}}, String}()
    buf = IOBuffer()
    first_residue = first(residue_list)
    key = (model=first_residue.id.model, chain=first_residue.id.chain)
    warnings = Dict{String, Set{String}}()
    for res in residue_list
        if !_is(res.id.model, model) || !_is(res.id.chain, chain) || !_is(res.id.group, group) || res.id.name == "HOH"
            continue
        end
        current_key = (model=res.id.model, chain=res.id.chain)
        if current_key != key
            _add_sequence!(buf, chains, key)
            key = current_key
        end
        if haskey(THREE2ONE, res.id.name)
            write(buf, THREE2ONE[res.id.name])
        else
            if haskey(warnings, res.id.name)
                push!(warnings[res.id.name], res.id.number)
            else
                warnings[res.id.name] = Set([res.id.number])
            end
        end
    end
    _add_sequence!(buf, chains, key)
    for (name, numbers) in warnings
        @warn("Residue name $(name) (residue numbers: $(join(numbers, ", "))) is not found in MIToS.Utils.THREE2ONE")
    end
    chains
end
