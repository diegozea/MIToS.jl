# Hobohm I
# ========

"""
Fill `cluster` and `clustersize` vectors. They are assumed to be empty (only
zeroes) and their length must be equal to the number of elements to cluster.
`within_cluster` is a predicate that takes two items and a `threshold` and
returns `true` if they should belong to the same cluster. `threshold` is passed
as the last argument to `within_cluster`. The number of elements is stored in
`n_items`.
"""
function _fill_hobohmI!(
    scan_candidates!::Function,
    state,
    within_cluster::Function,
    cluster::Vector{Int},
    clustersize::Vector{Int},
    items::AbstractVector,
    threshold,
)
    cluster_id = 0
    n_items = length(items)
    @inbounds for i = 1:(n_items-1)
        if cluster[i] == 0
            cluster_id += 1
            cluster[i] = cluster_id
            clustersize[cluster_id] += 1
            ref_item = items[i]
            clustersize[cluster_id] += scan_candidates!(
                state,
                within_cluster,
                cluster,
                items,
                ref_item,
                threshold,
                cluster_id,
                i + 1,
                n_items,
            )
        end
    end
    @inbounds if cluster[n_items] == 0
        cluster_id += 1
        cluster[n_items] = cluster_id
        clustersize[cluster_id] += 1
    end
    resize!(clustersize, cluster_id)
end

function _scan_hobohmI_serial!(
    ::Nothing,
    within_cluster::Function,
    cluster::Vector{Int},
    items::AbstractVector,
    ref_item,
    threshold,
    cluster_id::Int,
    first_candidate::Int,
    last_candidate::Int,
)
    added = 0
    @inbounds for j = first_candidate:last_candidate
        if cluster[j] == 0 && within_cluster(ref_item, items[j], threshold)
            cluster[j] = cluster_id
            added += 1
        end
    end
    added
end

function _scan_hobohmI_threaded!(
    matches::Vector{Bool},
    within_cluster::Function,
    cluster::Vector{Int},
    items::AbstractVector,
    ref_item,
    threshold,
    cluster_id::Int,
    first_candidate::Int,
    last_candidate::Int,
)
    Threads.@threads for j = first_candidate:last_candidate
        @inbounds matches[j] =
            cluster[j] == 0 && within_cluster(ref_item, items[j], threshold)
    end
    added = 0
    @inbounds for j = first_candidate:last_candidate
        if matches[j]
            cluster[j] = cluster_id
            added += 1
        end
    end
    added
end

function _fill_hobohmI!(
    within_cluster::Function,
    cluster::Vector{Int},
    clustersize::Vector{Int},
    items::AbstractVector,
    threshold;
    threads::Bool = true,
)
    use_threads = threads && Threads.nthreads() > 1
    scan_candidates! = ifelse(use_threads, _scan_hobohmI_threaded!, _scan_hobohmI_serial!)
    state = use_threads ? Vector{Bool}(undef, length(items)) : nothing
    _fill_hobohmI!(
        scan_candidates!,
        state,
        within_cluster,
        cluster,
        clustersize,
        items,
        threshold,
    )
end

"""
Calculates the weight of each sequence in a cluster. The weight is equal to one divided
by the number of sequences in the cluster.
"""
function _get_sequence_weight(clustersize, cluster)
    n_items = length(cluster)
    sequence_weight = Array{Float64}(undef, n_items)
    for i = 1:n_items
        @inbounds sequence_weight[i] = 1.0 / clustersize[cluster[i]]
    end
    Weights(sequence_weight, Float64(length(clustersize)))
end

function _hobohmI(within_cluster::Function, items::AbstractVector, threshold; threads::Bool)
    n = length(items)
    cluster = zeros(Int, n)
    clustersize = zeros(Int, n)
    _fill_hobohmI!(
        within_cluster,
        cluster,
        clustersize,
        items,
        threshold;
        threads = threads,
    )
    Clusters(clustersize, cluster, _get_sequence_weight(clustersize, cluster))
end

"""
`hobohmI(within_cluster, items, threshold; threads=false)`

Cluster `items` using the Hobohm I algorithm from Hobohm et al. `within_cluster`
is a predicate that receives two elements and `threshold` and returns `true` when
they should be clustered together. User-provided predicates remain serial by
default; set `threads` to `true` to thread the inner scan over candidate items
when worker threads are available.

# References

  - [Hobohm, Uwe, et al. "Selection of representative protein data sets."
    Protein Science 1.3 (1992): 409-417.](@cite 10.1002/pro.5560010313)
"""
function hobohmI(
    within_cluster::Function,
    items::AbstractVector,
    threshold;
    threads::Bool = false,
)
    _hobohmI(within_cluster, items, threshold; threads = threads)
end

function hobohmI(
    ::typeof(percentidentity),
    items::AbstractVector,
    threshold;
    threads::Bool = false,
)
    _hobohmI(percentidentity, items, threshold; threads = threads)
end

"""
`hobohmI(within_cluster, msa, threshold; threads=false)`

This method allows clustering the aligned sequences in `msa` using the
`within_cluster` predicate. It converts the alignment into a vector of
residue sequences and forwards the call to the general method. User-provided
predicates remain serial by default.
"""
function hobohmI(
    within_cluster::Function,
    msa::AbstractMatrix{Residue},
    threshold;
    threads::Bool = false,
)
    aln = getresiduesequences(msa)
    hobohmI(within_cluster, aln, threshold; threads = threads)
end

function hobohmI(
    ::typeof(percentidentity),
    msa::AbstractMatrix{Residue},
    threshold;
    threads::Bool = false,
)
    aln = getresiduesequences(msa)
    hobohmI(percentidentity, aln, threshold; threads = threads)
end

"""
`hobohmI(msa, threshold; threads=false)`

This method allows to cluster the sequences contained in `msa` using
`percentidentity` as the clustering predicate. When `threads` is `true`, the
inner scan is threaded if worker threads are available.
"""
function hobohmI(msa::AbstractMatrix{Residue}, threshold; threads::Bool = false)
    hobohmI(percentidentity, msa, threshold; threads = threads)
end
