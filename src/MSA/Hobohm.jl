# Hobohm I
# ========

const _HOBOHM_THREADS = """The `threads` keyword argument (default: `false`) controls 
whether the inner scan over candidate cluster members runs in parallel when worker threads 
are available."""

"""
Fill `cluster` with the Hobohm I assignments and return the number of clusters.
`cluster` is assumed to be empty (only zeroes) and its length must be equal to
the number of elements to cluster. `within_cluster` is a predicate that takes
two items and a `threshold` and returns `true` if they should belong to the
same cluster. `threshold` is passed as the last argument to `within_cluster`.
The number of elements is stored in `n_items`.
"""
function _fill_hobohmI!(
    scan_function::Function,
    within_cluster::Function,
    cluster::Vector{Int},
    items::AbstractVector,
    threshold,
)
    cluster_id = 0
    n_items = length(items)
    @inbounds for i = 1:(n_items-1)
        if cluster[i] == 0
            cluster_id += 1
            cluster[i] = cluster_id
            ref_item = items[i]
            scan_function(
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
    end
    cluster_id
end

function _scan_hobohmI_serial!(
    within_cluster::Function,
    cluster::Vector{Int},
    items::AbstractVector,
    ref_item,
    threshold,
    cluster_id::Int,
    first_candidate::Int,
    last_candidate::Int,
)
    @inbounds for j = first_candidate:last_candidate
        if cluster[j] == 0 && within_cluster(ref_item, items[j], threshold)
            cluster[j] = cluster_id
        end
    end
end

function _scan_hobohmI_threaded!(
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
        @inbounds if cluster[j] == 0 && within_cluster(ref_item, items[j], threshold)
            cluster[j] = cluster_id
        end
    end
end

function _fill_hobohmI!(
    within_cluster::Function,
    cluster::Vector{Int},
    items::AbstractVector,
    threshold;
    threads::Bool = true,
)
    use_threads = threads && Threads.nthreads() > 1
    scan_function = ifelse(use_threads, _scan_hobohmI_threaded!, _scan_hobohmI_serial!)
    _fill_hobohmI!(scan_function, within_cluster, cluster, items, threshold)
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
    nclusters = _fill_hobohmI!(within_cluster, cluster, items, threshold; threads = threads)
    resize!(clustersize, nclusters)
    @inbounds for i = 1:n
        clustersize[cluster[i]] += 1
    end
    Clusters(clustersize, cluster, _get_sequence_weight(clustersize, cluster))
end

"""
`hobohmI(within_cluster, items, threshold; threads=false)`

Cluster `items` using the Hobohm I algorithm from Hobohm et al. `within_cluster`
is a predicate that receives two elements and `threshold` and returns `true` when
they should be clustered together. $_HOBOHM_THREADS

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
residue sequences and forwards the call to the general method. $_HOBOHM_THREADS
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
`percentidentity` as the clustering predicate. $_HOBOHM_THREADS
"""
function hobohmI(msa::AbstractMatrix{Residue}, threshold; threads::Bool = false)
    hobohmI(percentidentity, msa, threshold; threads = threads)
end
