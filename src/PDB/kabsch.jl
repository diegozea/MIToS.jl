"""
    kabsch!(rotatedX::Matrix{Float64}, X::Matrix{Float64}, Y::Matrix{Float64})

Input: Two sets of points X,Y as NxD matrices, where D is dimension, N number of points, and a preallocated NxD matrix rotatedX.
Assumes that the centroids of X and Y are at the origin of coordinates.
Output: Rotates X so that RMSD(X,Y) is minimized. The rotated X is written into rotatedX.
"""
function kabsch!(rotatedX::Matrix{Float64}, X::Matrix{Float64}, Y::Matrix{Float64})
  @assert size(rotatedX) == size(X) == size(Y)
  A::Matrix{Float64} = X' * Y
  χ = eye(A)
  χ[end,end] = sign(det(A))
	u::Matrix{Float64}, σ::Vector{Float64}, v::Matrix{Float64} = svd(A)
  A_mul_B!(rotatedX, X, u * χ * v')
end

"""
    center!(X::Matrix{Float64})

Input: A set of points X as an NxD matrix (N: number of points, D: dimension)
Translates X in place so that its centroid is at the origin of coordinates
"""
function center!(X::Matrix{Float64})
  for i = 1:size(X)[2]
    X[:,i] -= mean(X[:,i])
  end
end

"""
    rmsd(A::Matrix{Float64}, B::Matrix{Float64})

Return RMSD between two sets of points A and B, given as NxD matrices (N: number of points, D: dimension)
"""
function rmsd(A::Matrix{Float64}, B::Matrix{Float64})
  @assert size(A) == size(B)
  N, D = size(A)
	s::Float64 = 0.0
	for i=1:N, j=1:D
		s += (A[i,j]::Float64 - B[i,j]::Float64)^2
	end
	return sqrt(s / N)
end
