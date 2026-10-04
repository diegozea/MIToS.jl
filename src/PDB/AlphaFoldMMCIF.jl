# A coordinate-derived template is deliberately separate from generic mmCIF export.
# PDBResidue does not retain the original polymer sequence or experimental metadata.

const _ALPHAFOLD_STANDARD_RESIDUES = (
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
)

# Peptide values in the wwPDB _chem_comp.type controlled vocabulary (case insensitive).
const _ALPHAFOLD_PEPTIDE_TYPES = (
    "peptide linking",
    "peptide-like",
    "l-peptide linking",
    "d-peptide linking",
    "l-peptide cooh carboxy terminus",
    "d-peptide cooh carboxy terminus",
    "l-peptide nh3 amino terminus",
    "d-peptide nh3 amino terminus",
    "l-beta-peptide, c-gamma linking",
    "d-beta-peptide, c-gamma linking",
    "l-gamma-peptide, c-delta linking",
    "d-gamma-peptide, c-delta linking",
)

function _alphafold_component_type(name::String, chem_comp_types::AbstractDict)
    value = if haskey(chem_comp_types, name)
        chem_comp_types[name]
    elseif name in _ALPHAFOLD_STANDARD_RESIDUES
        name == "GLY" ? "peptide linking" : "L-peptide linking"
    else
        throw(
            ArgumentError(
                "Residue $name needs an explicit peptide type in chem_comp_types. " *
                "Select only protein polymer residues; exclude waters, ligands and nucleic acids.",
            ),
        )
    end
    if !(value isa AbstractString) || !(lowercase(value) in _ALPHAFOLD_PEPTIDE_TYPES)
        throw(ArgumentError("Invalid peptide chem_comp type for $name: $value"))
    end
    String(value)
end

function _alphafold_residue_number(number::String)
    matched = match(r"^([+-]?[0-9]+)([A-Za-z]?)\z", number)
    if matched === nothing
        throw(
            ArgumentError(
                "Invalid author residue number $(repr(number)); expected an integer " *
                "with an optional single-letter insertion code.",
            ),
        )
    end
    num = tryparse(Int, matched[1])
    num === nothing &&
        throw(ArgumentError("Author residue number is out of range: $number"))
    (num, String(matched[2]))
end

"""
    alphafold_mmcifdict(residues; entry_id, model=nothing, exptl_method="?",
                       release_date=nothing, chem_comp_types=Dict(), renumber_auth=false)

Return a `BioStructures.MMCIFDict` for a coordinate-derived AlphaFold protein template.
Write it with `BioStructures.writemmcif(path, dict)`. Generic [`MMCIFFile`](@ref) exports
are unaffected, and the input residues and atoms are not modified.

Select only residues belonging to protein polymers. By passing them to this helper, the
caller declares their polymer membership. Waters, ligands and nucleic acids must be
excluded explicitly. Standard amino acids have known component types (glycine is
`"peptide linking"`); every other residue, including modified amino acids and `UNK`,
requires an explicit wwPDB peptide type in `chem_comp_types`, e.g.
`Dict("MSE" => "L-peptide linking")`. Names and `ATOM`/`HETATM` records are preserved.

Exactly one model must be selected. For multi-model input, supply its string identifier
with `model`. The selected model is written as model `1`, which AlphaFold expects.
Each chain gets a separate entity. Chains follow first appearance, and residues retain
their input order within each chain, even when author numbering is nonmonotonic.
`_entity_poly_seq` contains **only those observed residues**; missing residues and the
full experimental sequence cannot be recovered. New `label_seq_id`s run from `1` per
chain, replacing any stored PDBe numbering. Both label and author chain, component and
atom identifiers use the values stored in `PDBResidue`, which has only one namespace.

Author residue numbers must be integers with an optional single-letter insertion code.
Their numeric values and insertion codes are preserved by default. Duplicate author
positions, including microheterogeneity, are rejected. For ColabFold versions that skip
insertions, `renumber_auth=true` explicitly replaces author numbers with the new label
numbers and removes insertion codes **in the output only**. Empty chains, empty residues
and empty selections are rejected.

`entry_id` is a caller-chosen local identifier (letters, digits, `_`, `.`, `-`, starting
with a letter or digit). No experimental method, resolution, release date, missing
sequence or atom element is inferred. The method defaults to mmCIF unknown (`"?"`).
Supply a known `exptl_method` and, if available, `release_date::Dates.Date` from the
source structure. A supplied date is written as a single revision-history entry;
otherwise that category is omitted. AlphaFold date filtering may require it; ColabFold
may insert its own fallback date. This is a minimal template, not a wwPDB deposition or
a guarantee of acceptance by every downstream template pipeline.
"""
function alphafold_mmcifdict(
    residues::AbstractVector{PDBResidue};
    entry_id::AbstractString,
    model::Union{Nothing,AbstractString} = nothing,
    exptl_method::AbstractString = "?",
    release_date::Union{Nothing,Date} = nothing,
    chem_comp_types::AbstractDict = Dict{String,String}(),
    renumber_auth::Bool = false,
)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_.-]*\z", entry_id) ||
        throw(ArgumentError("entry_id must be a nonempty local mmCIF identifier"))
    isempty(strip(exptl_method)) &&
        throw(ArgumentError("Use '?' for an unknown experimental method"))

    selected =
        model === nothing ? collect(residues) : filter(r -> r.id.model == model, residues)
    isempty(selected) && throw(ArgumentError("No residues selected for the template"))
    models = unique(r.id.model for r in selected)
    length(models) == 1 ||
        throw(ArgumentError("Select exactly one template model with the model keyword"))

    chains = OrderedDict{String,Vector{PDBResidue}}()
    component_types = OrderedDict{String,String}()
    for res in selected
        chain = res.id.chain
        if isempty(chain) || chain in (".", "?") || any(isspace, chain)
            throw(
                ArgumentError("Assign a nonempty chain identifier before template export"),
            )
        end
        isempty(res.atoms) && throw(ArgumentError("Template residues must contain atoms"))
        res.id.group in ("ATOM", "HETATM") ||
            throw(ArgumentError("Invalid atom record group: $(res.id.group)"))
        if isempty(res.id.name) || res.id.name in (".", "?") || any(isspace, res.id.name)
            throw(ArgumentError("Template residues must have a component identifier"))
        end
        component_types[res.id.name] =
            _alphafold_component_type(res.id.name, chem_comp_types)
        push!(get!(Vector{PDBResidue}, chains, chain), res)
    end

    template_residues = PDBResidue[]
    atom_entities = String[]
    seq_entities = String[]
    seq_numbers = String[]
    seq_monomers = String[]
    for (entity_index, (chain, chain_residues)) in enumerate(chains)
        entity = string(entity_index)
        positions = Set{Tuple{Int,String}}()
        for (index, res) in enumerate(chain_residues)
            number, inscode = _alphafold_residue_number(res.id.number)
            position = (number, inscode)
            position in positions && throw(
                ArgumentError(
                    "Duplicate author position $(res.id.number) in chain $chain; " *
                    "resolve duplicate residues or microheterogeneity before export",
                ),
            )
            push!(positions, position)
            seq_number = string(index)
            auth_number = renumber_auth ? seq_number : string(number, inscode)
            id = PDBResidueIdentifier(
                seq_number,
                auth_number,
                res.id.name,
                res.id.group,
                "1",
                chain,
            )
            push!(template_residues, PDBResidue(id, res.atoms))
            append!(atom_entities, fill(entity, length(res.atoms)))
            push!(seq_entities, entity)
            push!(seq_numbers, seq_number)
            push!(seq_monomers, res.id.name)
        end
    end

    mmcif = _pdbresidues_to_mmcifdict(template_residues; label = true)
    for field in ("asym_id", "comp_id", "atom_id")
        mmcif["_atom_site.auth_$field"] = copy(mmcif["_atom_site.label_$field"])
    end
    mmcif["_atom_site.label_entity_id"] = atom_entities
    mmcif["data_"] = [String(entry_id)]
    mmcif["_entry.id"] = [String(entry_id)]
    mmcif["_exptl.entry_id"] = [String(entry_id)]
    mmcif["_exptl.method"] = [String(exptl_method)]
    entities = string.(1:length(chains))
    mmcif["_entity.id"] = entities
    mmcif["_entity.type"] = fill("polymer", length(chains))
    mmcif["_struct_asym.id"] = collect(keys(chains))
    mmcif["_struct_asym.entity_id"] = copy(entities)
    mmcif["_entity_poly_seq.entity_id"] = seq_entities
    mmcif["_entity_poly_seq.num"] = seq_numbers
    mmcif["_entity_poly_seq.mon_id"] = seq_monomers
    mmcif["_entity_poly_seq.hetero"] = fill("n", length(seq_numbers))
    mmcif["_chem_comp.id"] = collect(keys(component_types))
    mmcif["_chem_comp.type"] = collect(values(component_types))
    if release_date !== nothing
        mmcif["_pdbx_audit_revision_history.ordinal"] = ["1"]
        mmcif["_pdbx_audit_revision_history.revision_date"] = [string(release_date)]
    end
    mmcif
end
