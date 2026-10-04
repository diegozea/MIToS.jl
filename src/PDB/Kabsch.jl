# Jorge Fernández de Cossío Díaz ( @cossio ) wrote the kabsch, rmsd and center! functions.

"""
`kabsch(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64})`

This function takes two sets of points, `A` (refrence) and `B` as NxD matrices, where D
is the dimension and N is the number of points.
Assumes that the centroids of `A` and `B` are at the origin of coordinates.
You can call `center!` on each matrix before calling `kabsch` to center the matrices
in the `(0.0, 0.0, 0.0)`.
Rotates `B` so that `rmsd(A,B)` is minimized.
Returns the rotation matrix. You should do `B * RotationMatrix` to get the rotated B.
"""
function kabsch(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64})
    @assert size(A) == size(B)
    M::AbstractMatrix{Float64} = B' * A
    χ = Matrix{Float64}(I, size(M, 1), size(M, 2))
    u::AbstractMatrix{Float64}, σ::Vector{Float64}, v::AbstractMatrix{Float64} = svd(M)
    # det(M) can be zero for planar/collinear fragments. The orthogonal factors
    # still define a proper rotation; a zero correction would collapse one axis.
    χ[end, end] = det(u * v') < 0.0 ? -1.0 : 1.0
    return u * χ * v'
end

"""
`center!(A::AbstractMatrix{Float64})`

Takes a set of points `A` as an NxD matrix (N: number of points, D: dimension).
Translates `A` in place so that its centroid is at the origin of coordinates
"""
function center!(A::AbstractMatrix{Float64})
    for i = 1:size(A)[2]
        A[:, i] .-= mean(A[:, i])
    end
end

"""
`rmsd(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64})`

Return RMSD between two sets of points `A` and `B`, given as NxD matrices
(N: number of points, D: dimension).
"""
function rmsd(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64})
    @assert size(A) == size(B)
    N, D = size(A)
    s::Float64 = 0.0
    for i = 1:N, j = 1:D
        s += (B[i, j]::Float64 - A[i, j]::Float64)^2
    end
    return sqrt(s / N)
end

"""
`rmsd(A::AbstractVector{PDBResidue}, B::AbstractVector{PDBResidue}; superimposed::Bool=false)`

Returns the Cα RMSD value between two PDB structures: `A` and `B`.
If the structures are already superimposed between them,
use `superimposed=true` to avoid a new superimposition (`superimposed` is `false` by default).
"""
function rmsd(
    A::AbstractVector{PDBResidue},
    B::AbstractVector{PDBResidue};
    superimposed::Bool = false,
)
    if superimposed
        rmsd(CAmatrix(A), CAmatrix(B))
    else
        superimpose(A, B)[end]::Float64
    end
end

# Kabsch for Vector{PDBResidues}
# ==============================

"""
Returns the Cα with best occupancy in the `PDBResidue`.
If the `PDBResidue` has no Cα, `missing` is returned.
"""
function getCA(res::PDBResidue)
    if length(res) == 0
        @warn """There are no atoms in residue
        $(res.id)"""
        return missing
    end
    CAs = findatoms(res, "CA")
    if length(CAs) == 0
        @warn """There is no alpha carbon in residue 
        $(res.id)"""
        missing
    else
        CAindex = selectbestoccupancy(res, CAs)
        res.atoms[CAindex]
    end
end

"""
Returns a matrix with the x, y and z coordinates of the Cα with best occupancy for each
`PDBResidue` of the ATOM group. If a residue doesn't have a Cα, its Cα coordinates are NaNs.
"""
function CAmatrix(residues::AbstractVector{PDBResidue})
    len = length(residues)
    CAlist = Array{Float64}(undef, 3 * len)
    j = 0
    r = 0
    @inbounds for i = 1:len
        res = residues[i]
        if (res.id.group == "ATOM") && (length(res) > 0)
            r += 1
            CAs = findatoms(res, "CA")
            if length(CAs) != 0
                CAindex = selectbestoccupancy(res, CAs)
                coord = res.atoms[CAindex].coordinates
                CAlist[j+=1] = coord.x
                CAlist[j+=1] = coord.y
                CAlist[j+=1] = coord.z
            else
                CAlist[j+=1] = NaN
                CAlist[j+=1] = NaN
                CAlist[j+=1] = NaN
            end
        end
    end
    reshape(resize!(CAlist, j), (3, r))'
end

"""
Returns a matrix with the x, y, z coordinates of each atom in each `PDBResidue`
"""
function coordinatesmatrix(res::PDBResidue)
    atoms = res.atoms
    len = length(atoms)
    mat = Array{Float64}(undef, 3, len)
    for i = 1:len
        coord = atoms[i].coordinates
        mat[1, i] = coord.x
        mat[2, i] = coord.y
        mat[3, i] = coord.z
    end
    mat'
end

function coordinatesmatrix(residues::AbstractVector{PDBResidue})
    reduce(vcat, map(coordinatesmatrix, residues))
end

"""
Returns a `Matrix{Float64}` with the centered coordinates of all the atoms in `residues`.
An optional positional argument `CA` (default: `true`) defines if only Cα carbons should
be used to center the matrix.
"""
function centeredcoordinates(residues::AbstractVector{PDBResidue}, CA::Bool = true)
    coordinates = PDB.coordinatesmatrix(residues)
    meancoord = CA ? mean(PDB.CAmatrix(residues), dims = 1) : mean(coordinates, dims = 1)
    coordinates .- meancoord
end

"""
Returns a new `Vector{PDBResidue}` with the `PDBResidue`s having centered coordinates.
An optional positional argument `CA` (default: `true`) defines if only Cα carbons should
be used to center the matrix.
"""
function centeredresidues(residues::AbstractVector{PDBResidue}, CA::Bool = true)
    coordinates = centeredcoordinates(residues, CA)
    change_coordinates(residues, coordinates)
end

"""
`change_coordinates(atom::PDBAtom, coordinates::Coordinates)`

Returns a new `PDBAtom` but with a new `coordinates`
"""
function change_coordinates(atom::PDBAtom, coordinates::Coordinates)
    PDBAtom(
        coordinates,
        identity(atom.atom),
        identity(atom.element),
        copy(atom.occupancy),
        identity(atom.B),
        identity(atom.alt_id),
        identity(atom.charge),
    )
end

"""
`change_coordinates(residue::PDBResidue, coordinates::AbstractMatrix{Float64}, offset::Int=1)`

Returns a new `PDBResidues` with (x,y,z) from a coordinates `AbstractMatrix{Float64}`
You can give an `offset` indicating in which matrix row starts the (x,y,z) coordinates
of the residue.
"""
function change_coordinates(
    residue::PDBResidue,
    coordinates::AbstractMatrix{Float64},
    offset::Int = 1,
)
    centeredatoms = map(residue.atoms) do atom
        atoms = change_coordinates(atom, Coordinates(vec(coordinates[offset, :])))
        offset += 1
        return atoms
    end
    PDBResidue(residue.id, centeredatoms)
end

"""
`change_coordinates(residues::AbstractVector{PDBResidue}, coordinates::AbstractMatrix{Float64})`

Returns a new `Vector{PDBResidues}` with (x,y,z) from a coordinates `Matrix{Float64}`
"""
function change_coordinates(
    residues::AbstractVector{PDBResidue},
    coordinates::AbstractMatrix{Float64},
)
    nres = length(residues)
    updated = Array{PDBResidue}(undef, nres)
    j = 1
    for i = 1:nres
        residue = residues[i]
        updated[i] = change_coordinates(residue, coordinates, j)
        j += length(residue)
    end
    updated
end

_iscentered(x::Float64, y::Float64, z::Float64) =
    (abs(x) <= 1e-13) && (abs(y) <= 1e-13) && (abs(z) <= 1e-13)

_iscentered(meanCα::AbstractVector{Float64}) = _iscentered(meanCα[1], meanCα[2], meanCα[3])

_iscentered(CA::AbstractMatrix{Float64}) = _iscentered(vec(mean(CA, dims = 1)))


"""
Return the matching CA matrices after deleting the rows/residues where the CA
is missing in at least one structure.
"""
function _get_matched_Cαs(
    A::AbstractVector{PDBResidue},
    B::AbstractVector{PDBResidue},
    ::Nothing,
)
    length_A = length(A)
    @assert length_A == length(B) "PDBResidue vectors should have the same length."
    ACα = PDB.CAmatrix(A)
    BCα = PDB.CAmatrix(B)
    without_Cα = isnan.(ACα[:, 1]) .| isnan.(BCα[:, 1])
    if any(without_Cα)
        n_without_ca = sum(without_Cα)
        @assert length_A - n_without_ca != 0 "There are not alpha-carbons to align."
        @warn string(
            "Using ",
            length_A - n_without_ca,
            " residues for RMSD calculation because there are ",
            n_without_ca,
            " residues without CA: ",
            findall(without_Cα),
        )
        with_Cα = .!without_Cα
        return ACα[with_Cα, :], BCα[with_Cα, :]
    end
    ACα, BCα
end

function _get_matched_Cαs(
    A::AbstractVector{PDBResidue},
    B::AbstractVector{PDBResidue},
    matches,
)
    if Base.IteratorSize(typeof(matches)) == Base.SizeUnknown()
        Asel, Bsel = PDBResidue[], PDBResidue[]
        for (i, j) in matches
            push!(Asel, A[i])
            push!(Bsel, B[j])
        end
        return _get_matched_Cαs(Asel, Bsel, nothing)
    end
    Asel = Vector{PDBResidue}(undef, length(matches))
    Bsel = similar(Asel)
    for (k, (i, j)) in enumerate(matches)
        Asel[k] = A[i]
        Bsel[k] = B[j]
    end
    return _get_matched_Cαs(Asel, Bsel, nothing)
end

"""
    Asuper, Bsuper, RMSD = superimpose(A, B, matches=nothing)

This function takes `A::AbstractVector{PDBResidue}` (reference) and
`B::AbstractVector{PDBResidue}`. Translates `A` and `B` to the origin of coordinates,
and rotates `B` so that `rmsd(A,B)` is minimized with the Kabsch algorithm
(using only their α carbons).
Returns the rotated and translated versions of `A` and `B`, and the RMSD value.

Optionally provide `matches` which iterates over matched index pairs in `A` and `B`,
e.g., `matches = [(3, 5), (4, 6), ...]`. The alignment will be constructed
using just the matching residues.
"""
function superimpose(
    A::AbstractVector{PDBResidue},
    B::AbstractVector{PDBResidue},
    matches = nothing,
)
    ACα, BCα = _get_matched_Cαs(A, B, matches)
    Bxyz = PDB.coordinatesmatrix(B)
    meanACα = vec(mean(ACα, dims = 1))
    meanBCα = vec(mean(BCα, dims = 1))
    if !_iscentered(meanBCα)
        @inbounds BCα[:, 1] .-= meanBCα[1]
        @inbounds BCα[:, 2] .-= meanBCα[2]
        @inbounds BCα[:, 3] .-= meanBCα[3]
        @inbounds Bxyz[:, 1] .-= meanBCα[1]
        @inbounds Bxyz[:, 2] .-= meanBCα[2]
        @inbounds Bxyz[:, 3] .-= meanBCα[3]
    end
    if !_iscentered(meanACα)
        @inbounds ACα[:, 1] .-= meanACα[1]
        @inbounds ACα[:, 2] .-= meanACα[2]
        @inbounds ACα[:, 3] .-= meanACα[3]
        Axyz = PDB.coordinatesmatrix(A)
        @inbounds Axyz[:, 1] .-= meanACα[1]
        @inbounds Axyz[:, 2] .-= meanACα[2]
        @inbounds Axyz[:, 3] .-= meanACα[3]
        RotationMatrix = PDB.kabsch(ACα, BCα)
        return (
            change_coordinates(A, Axyz),
            change_coordinates(B, Bxyz * RotationMatrix),
            PDB.rmsd(ACα, BCα * RotationMatrix),
        )
    else
        RotationMatrix = PDB.kabsch(ACα, BCα)
        return (
            A,
            change_coordinates(B, Bxyz * RotationMatrix),
            PDB.rmsd(ACα, BCα * RotationMatrix),
        )
    end
end

# Structure similarity metrics
# ----------------------------

_default_matches(A::AbstractVector{PDBResidue}, B::AbstractVector{PDBResidue}) =
    ((i, i) for i = 1:min(length(A), length(B)))

@inline _has_ca(res::PDBResidue) = !isempty(findatoms(res, "CA"))

function _valid_matches(
    A::AbstractVector{PDBResidue},
    B::AbstractVector{PDBResidue},
    matches,
)
    Base.require_one_based_indexing(A, B)
    valid = Tuple{Int,Int}[]
    used_A, used_B = Set{Int}(), Set{Int}()
    for (i, j) in matches
        (i isa Integer && j isa Integer) ||
            throw(ArgumentError("matches must contain integer index pairs"))
        checkbounds(A, i)
        checkbounds(B, j)
        (i in used_A || j in used_B) && throw(
            ArgumentError("matches must be one-to-one; each index may occur only once"),
        )
        push!(used_A, i)
        push!(used_B, j)
        if _has_ca(A[i]) && _has_ca(B[j])
            push!(valid, (i, j))
        end
    end
    return valid
end

function _similarity_coordinates(A, B, matches, Ltarget, max_iterations)
    Ltarget >= 0 || throw(ArgumentError("Ltarget must be nonnegative"))
    max_iterations >= 0 || throw(ArgumentError("max_iterations must be nonnegative"))
    matches_iter = isnothing(matches) ? _default_matches(A, B) : matches
    pairs = _valid_matches(A, B, matches_iter)
    n = length(pairs)
    Ltarget >= n || throw(ArgumentError("Ltarget must be at least the number of Cα pairs"))
    ACα, BCα = Matrix{Float64}(undef, n, 3), Matrix{Float64}(undef, n, 3)
    for (k, (i, j)) in enumerate(pairs)
        a, b = getCA(A[i]).coordinates, getCA(B[j]).coordinates
        (all(isfinite, a) && all(isfinite, b)) ||
            throw(ArgumentError("matched Cα coordinates must be finite"))
        ACα[k, :] = a
        BCα[k, :] = b
    end
    return ACα, BCα
end

# Fit only the selected Cα pairs, then measure every pair under that same rigid fit.
# Keeping coordinates here avoids repeatedly copying all atoms in both structures.
function _ca_distances(A::Matrix{Float64}, B::Matrix{Float64}, subset)
    Asel, Bsel = A[subset, :], B[subset, :]
    mean_A, mean_B = mean(Asel, dims = 1), mean(Bsel, dims = 1)
    rotation = kabsch(Asel .- mean_A, Bsel .- mean_B)
    delta = (B .- mean_B) * rotation .- (A .- mean_A)
    return vec(sqrt.(sum(abs2, delta, dims = 2)))
end

function _fragment_windows(n::Int, window_sizes)
    windows = collect(window_sizes)
    all(w -> w isa Integer && w >= 3, windows) ||
        throw(ArgumentError("window_sizes must contain integers of at least 3"))
    fragments = UnitRange{Int}[1:n]
    for w in unique(windows)
        w >= n && continue
        # Visit every starting position: striding by half a window can miss a core.
        for start = 1:(n-w+1)
            push!(fragments, start:(start+Int(w)-1))
        end
    end
    return fragments
end

function _update_gdt_counts!(best_counts, distances, cutoffs)
    for k in eachindex(cutoffs)
        best_counts[k] = max(best_counts[k], count(d -> d <= cutoffs[k], distances))
    end
end

function _best_counts(A, B, fragments, cutoffs, max_iterations)
    best_counts = zeros(Int, length(cutoffs))
    for fragment in fragments
        initial_distances = _ca_distances(A, B, fragment)
        _update_gdt_counts!(best_counts, initial_distances, cutoffs)
        for cutoff in cutoffs
            subset = fragment
            distances = initial_distances
            for _ = 1:max_iterations
                next_subset = findall(d -> d <= cutoff, distances)
                (isempty(next_subset) || next_subset == subset) && break
                distances = _ca_distances(A, B, next_subset)
                # Retain the best count at every cutoff, not just the final fit.
                _update_gdt_counts!(best_counts, distances, cutoffs)
                subset = next_subset
            end
        end
    end
    return best_counts
end

"""
    gdt_per_cutoff(A, B; matches = nothing, cutoffs = (1.0, 2.0, 4.0, 8.0),
                  Ltarget = length(A), local_search = true,
                  window_sizes = (4, 8, 16, 32), max_iterations = 20)

Approximate the Global Distance Test (GDT) percentages between protein residues `A`
(reference) and `B` (model). Return a `Dict{Float64,Float64}` mapping each distance
cutoff in Å to `100 * N / Ltarget`, where `N` is the largest number of Cα pairs found
within that cutoff (inclusive). Each cutoff can use a different superposition.

`matches` iterates over one-to-one index pairs `(i, j)` in `A` and `B`. By default,
positions up to the shorter vector are paired; no sequence alignment or matching by
residue number is performed. Out-of-bounds or repeated indices are errors. Pairs
missing a Cα are skipped, without reducing the denominator. The highest-occupancy Cα
is used, including for modified amino acids stored as HETATM. Select the protein
residues/chain/model of interest before calling this function. Nonfinite matched
coordinates are errors. Inputs are not modified.

`Ltarget` is the reference length, including unpaired residues (default `length(A)`),
and must be at least the number of usable pairs. Thus incomplete models are penalized.
Pass the length of the intended evaluation domain explicitly when needed. No usable
pairs, including empty structures with `Ltarget = 0`, give zero scores.

With `local_search = true`, the whole alignment and every contiguous window in the
ordered list of usable pairs seed separate iterative fits for each cutoff. Each fit
repeatedly superposes the pairs within that cutoff, stopping when the subset is
unchanged/empty or after `max_iterations` refits. All intermediate scores are retained.
This is a fragment-search heuristic inspired by *Zemla*, not an exact LGA implementation
or a guaranteed global optimum. `local_search = false` evaluates just one global
least-squares fit. `max_iterations = 0` evaluates seeds without refinement.
`window_sizes` must contain integers ≥ 3; `cutoffs` must be finite and positive.

# References

  - [Zemla, Adam. "LGA: a method for finding 3D similarities in protein structures."
    Nucleic Acids Research 31.13 (2003): 3370–3374.](@cite 10.1093/nar/gkg571)
"""
function gdt_per_cutoff(
    A::AbstractVector{PDBResidue},
    B::AbstractVector{PDBResidue};
    matches = nothing,
    cutoffs = (1.0, 2.0, 4.0, 8.0),
    Ltarget::Integer = length(A),
    local_search::Bool = true,
    window_sizes = (4, 8, 16, 32),
    max_iterations::Integer = 20,
)
    cutoffs_vec = unique(Float64.(collect(cutoffs)))
    all(c -> isfinite(c) && c > 0, cutoffs_vec) ||
        throw(ArgumentError("cutoffs must be finite and positive"))
    ACα, BCα = _similarity_coordinates(A, B, matches, Ltarget, max_iterations)
    n = size(ACα, 1)
    fragments = _fragment_windows(n, window_sizes)
    if n == 0 || isempty(cutoffs_vec)
        return Dict(c => 0.0 for c in cutoffs_vec)
    end
    if !local_search
        fragments = [1:n]
    end
    best_counts =
        _best_counts(ACα, BCα, fragments, cutoffs_vec, local_search ? max_iterations : 0)
    return Dict(c => 100.0 * best_counts[k] / Ltarget for (k, c) in enumerate(cutoffs_vec))
end

"""
    gdt_ts(A, B; kwargs...)

Approximate GDT_TS on a 0–100 scale: the mean of independently maximized GDT
percentages at 1, 2, 4 and 8 Å. Keyword arguments (except `cutoffs`, which are fixed)
are forwarded to [`gdt_per_cutoff`](@ref), including reference-length normalization
and the fragment-search limitations described there.
"""
function gdt_ts(A::AbstractVector{PDBResidue}, B::AbstractVector{PDBResidue}; kwargs...)
    haskey(kwargs, :cutoffs) &&
        throw(ArgumentError("use gdt_per_cutoff for custom cutoffs"))
    scores = gdt_per_cutoff(A, B; cutoffs = (1.0, 2.0, 4.0, 8.0), kwargs...)
    return mean(values(scores))
end

"""
    gdt_ha(A, B; kwargs...)

Approximate GDT_HA on a 0–100 scale: the mean of independently maximized GDT
percentages at 0.5, 1, 2 and 4 Å. Keyword arguments (except `cutoffs`, which are fixed)
are forwarded to [`gdt_per_cutoff`](@ref), including reference-length normalization
and the fragment-search limitations described there.
"""
function gdt_ha(A::AbstractVector{PDBResidue}, B::AbstractVector{PDBResidue}; kwargs...)
    haskey(kwargs, :cutoffs) &&
        throw(ArgumentError("use gdt_per_cutoff for custom cutoffs"))
    scores = gdt_per_cutoff(A, B; cutoffs = (0.5, 1.0, 2.0, 4.0), kwargs...)
    return mean(values(scores))
end

# Convert before subtracting to also support unsigned reference lengths.
_tm_d0(L::Integer) = max(0.5, 1.24 * cbrt(Float64(L) - 15.0) - 1.8)

function _tm_score_from_distances(distances::AbstractVector{<:Real}, Ltarget::Integer)
    d0 = _tm_d0(Ltarget)
    s = 0.0
    @inbounds for d in distances
        s += 1.0 / (1.0 + (d / d0)^2)
    end
    return s / Ltarget
end

function _tm_subset(distances, cutoff)
    # At least three pairs (or all pairs for shorter inputs) determine a refit.
    # Expanding directly to the third-nearest pair also avoids unbounded loops.
    min_distance = partialsort(distances, min(3, length(distances)))
    threshold = max(cutoff, min_distance)
    return findall(d -> d <= threshold, distances)
end

"""
    tm_score(A, B; matches = nothing, Ltarget = length(A), local_search = true,
             window_sizes = (4, 8, 16, 32), max_iterations = 20)

Approximate the TM-score of protein model `B` against reference `A`, using the
*Zhang and Skolnick* sum `sum(1 / (1 + (dᵢ / d₀)^2)) / Ltarget` over all usable Cα
pairs. Here `d₀ = max(0.5, 1.24 * cbrt(Ltarget - 15) - 1.8)` Å. The reference length
sets both the denominator and distance scale, so swapping `A` and `B` can change the
score. Scores lie in [0, 1]. Missing/unpaired residues contribute zero.

Pairing, coordinate selection, normalization and empty-input handling are as in
[`gdt_per_cutoff`](@ref). In particular, `matches` is a supplied correspondence,
not an instruction to find a sequence or structural alignment.

With `local_search = true`, global and contiguous fragment fits seed iterative
superposition of nearby pairs. Following the TM-score search, the selection radius
is `clamp(d₀, 4.5, 8.0) - 1` Å initially and `clamp(d₀, 4.5, 8.0) + 1` Å for refits,
expanded when necessary to include at least three pairs. Every fit is scored over
all pairs, including those outside the selection radius, and the maximum is retained.
Refinement stops at an unchanged subset or `max_iterations` (default 20).
This heuristic is not guaranteed to reproduce the standalone TM-score program or
its optimum. `local_search = false` scores one global least-squares fit;
`max_iterations = 0` scores only the seeds.

# References

  - [Zhang, Yang, and Jeffrey Skolnick. "Scoring function for automated assessment of
    protein structure template quality." Proteins 57.4 (2004): 702–710.](@cite 10.1002/prot.20264)
"""
function tm_score(
    A::AbstractVector{PDBResidue},
    B::AbstractVector{PDBResidue};
    matches = nothing,
    Ltarget::Integer = length(A),
    local_search::Bool = true,
    window_sizes = (4, 8, 16, 32),
    max_iterations::Integer = 20,
)
    ACα, BCα = _similarity_coordinates(A, B, matches, Ltarget, max_iterations)
    n = size(ACα, 1)
    fragments = _fragment_windows(n, window_sizes)
    n == 0 && return 0.0
    if !local_search
        fragments = [1:n]
    end
    search_radius = clamp(_tm_d0(Ltarget), 4.5, 8.0)
    best_score = 0.0
    for fragment in fragments
        distances = _ca_distances(ACα, BCα, fragment)
        best_score = max(best_score, _tm_score_from_distances(distances, Ltarget))
        subset = fragment
        cutoff = search_radius - 1.0
        for _ = 1:(local_search ? max_iterations : 0)
            next_subset = _tm_subset(distances, cutoff)
            # The radius grows after the first selection, so even an unchanged
            # initial subset must be evaluated once with the larger radius.
            next_subset == subset && cutoff == search_radius + 1.0 && break
            distances = _ca_distances(ACα, BCα, next_subset)
            best_score = max(best_score, _tm_score_from_distances(distances, Ltarget))
            subset = next_subset
            cutoff = search_radius + 1.0
        end
    end
    return best_score
end

# RMSF: Root Mean-Square-average distance (Fluctuation)
# -----------------------------------------------------

"""
This looks for errors in the input to rmsf methods
"""
function _rmsf_test(vector)
    n = length(vector)
    @assert n >= 2 "You need at least two matrices/structures"
    sizes = unique(Tuple{Int,Int}[size(s) for s in vector])
    @assert length(sizes) == 1 "Matrices/Structures must have the same number of rows/atoms"
    @assert sizes[1][2] == 3 "Matrices should have 3 columns (x, y, z)"
end

"""
Calculates the average/mean position of each atom in a set of structure.
The function takes a vector (`AbstractVector`) of vectors (`AbstractVector{PDBResidue}`)
or matrices (`AbstractMatrix{Float64}`) as first argument. As second (optional) argument this
function can take an `AbstractVector{Float64}` of matrix/structure weights to return a
weighted mean. When a AbstractVector{PDBResidue} is used, if the keyword argument `calpha` is
`false` the RMSF is calculated for all the atoms. By default only alpha carbons are used
(default: `calpha=true`).
"""
function mean_coordinates(vec::AbstractVector{T}) where {T<:AbstractMatrix{Float64}}
    _rmsf_test(vec)
    n = length(vec)
    reduce(+, vec) ./ n
end

function mean_coordinates(
    vec::AbstractVector{T},
    matrixweights::AbstractVector{Float64},
) where {T<:AbstractMatrix{Float64}}
    _rmsf_test(vec)
    @assert length(vec) == length(matrixweights) "The number of matrix weights must be equal to the number of matrices."
    n = sum(matrixweights)
    reduce(+, (vec .* matrixweights)) ./ n
end

function mean_coordinates(
    vec::AbstractVector{T};
    calpha::Bool = true,
) where {T<:AbstractVector{PDBResidue}}
    mean_coordinates(map(calpha ? CAmatrix : coordinatesmatrix, vec))
end

function mean_coordinates(
    vec::AbstractVector{T},
    args...;
    calpha::Bool = true,
) where {T<:AbstractVector{PDBResidue}}
    mean_coordinates(map(calpha ? CAmatrix : coordinatesmatrix, vec), args...)
end

"""
Calculates the RMSF (Root Mean-Square-Fluctuation) between an atom and its average
position in a set of structures. The function takes a vector (`AbstractVector`) of
vectors (`AbstractVector{PDBResidue}`) or matrices (`AbstractMatrix{Float64}`) as first
argument. As second (optional) argument this function can take an `AbstractVector{Float64}`
of matrix/structure weights to return the root weighted mean-square-fluctuation around
the weighted mean structure. When a Vector{PDBResidue} is used, if the keyword argument
`calpha` is `false` the RMSF is calculated for all the atoms. By default only alpha
carbons are used (default: `calpha=true`).
"""
function rmsf(vector::AbstractVector{T}) where {T<:AbstractMatrix{Float64}}
    m = mean_coordinates(vector)
    # MIToS RMSF is calculated as in Eq. 6 from:
    # Kuzmanic, Antonija, and Bojan Zagrovic.
    # "Determination of ensemble-average pairwise root mean-square deviation from experimental B-factors."
    # Biophysical journal 98.5 (2010): 861-871.
    vec(sqrt.(mean(map(mat -> mapslices(x -> sum(abs2, x), mat .- m, dims = 2), vector))))
end

function rmsf(
    vector::AbstractVector{T},
    matrixweights::AbstractVector{Float64},
) where {T<:AbstractMatrix{Float64}}
    m = mean_coordinates(vector, matrixweights)
    d = map(mat -> mapslices(x -> sum(abs2, x), mat .- m, dims = 2), vector)
    vec(sqrt.(sum(d .* matrixweights) / sum(matrixweights)))
end

function rmsf(
    vec::AbstractVector{T};
    calpha::Bool = true,
) where {T<:AbstractVector{PDBResidue}}
    rmsf(map(calpha ? CAmatrix : coordinatesmatrix, vec))
end

function rmsf(
    vec::AbstractVector{T},
    matrixweights::AbstractVector{Float64};
    calpha::Bool = true,
) where {T<:AbstractVector{PDBResidue}}
    rmsf(map(calpha ? CAmatrix : coordinatesmatrix, vec), matrixweights)
end
