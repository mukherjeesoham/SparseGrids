#-----------------------------------------------------
# SparseGrids | Hierarchical Basis
# Soham 07/2020
#-----------------------------------------------------

using LinearAlgebra
export Level, HierarchicalBasis, NodalBasis
export basis, x, levels

#---------------------------------------------------------------
# Construct basis function and collocation points
#---------------------------------------------------------------

function basis(l::Int, j::CartesianIndex{D}, x::NTuple{D,T})::T where {T,D}
    """
    Evaluate basis functions for level (l) and position (j).
    range(l) := [0, ...]
    range(j) := [0, 2^l]
    """
    j = Tuple(j) .- 1
	@assert 0 <= l
    @assert all(0 .<= j .<= (2^l))

    h = 1/(2^l)

	if all(h.*(j .- 1) .<= x .<= h.*(j .+ 1)) && all(0 .<= x .<= 1)
        value = prod(1 .- abs.((x./h) .- j))
	else
		value = 0
	end
	return value
end

function x(l::Int, j::CartesianIndex{D})::NTuple{D,Float64} where {D}
    j = Tuple(j)
	@assert 0 <= l
    @assert all(1 .<= j .<= (2^l)+1)
    h = 1/(2^l)
    return h.*(j .- 1)
end

#---------------------------------------------------------------
# Construct types for arbitrary dimensions
# Assume we have equal number of levels and points in
# all directions.
#---------------------------------------------------------------

struct Level{T,D}
    values::Array{T,D}
end
    
struct HierarchicalBasis{T,D}
    levels::Array{Level{T,D},1}
end

struct NodalBasis{T,D}
    values::Array{T,D}
end

#---------------------------------------------------------------
# Construct indexing utilities
#---------------------------------------------------------------

function levels(HB::HierarchicalBasis{T,D})::Int where {T,D}
    return length(HB.levels) - 1
end

function levels(NB::NodalBasis{T,D})::Int where {T,D}
    return log2(size(NB.values,1) - 1)
end

function Base.getindex(HB::HierarchicalBasis{T,D}, l::Int, j::CartesianIndex{D})::T where {T,D}
    @assert 0  <= l  <= levels(HB)
    @assert all(1 .<= Tuple(j) .<= (2^l)+1)
    return HB.levels[l+1].values[j]
end

function Base.getindex(NB::NodalBasis{T}, j::CartesianIndex{D})::T where {T,D}
    @assert all(1 .<= Tuple(j) .<= (2^levels(NB))+1)
    return NB.values[j] 
end

function Base.getindex(NB::NodalBasis{T}, l::Int, j::CartesianIndex{D})::T where {T,D}
    @assert all(1 .<= Tuple(j) .<= (2^levels(NB))+1)
    return NB.values[(2^(levels(NB)-l))*j] 
end

function Base.setindex!(HB::HierarchicalBasis{T}, X::T, l::Int, j::CartesianIndex{D})::T where {T,D}
    @assert 0  <= l  <= levels(HB)
    @assert all(1 .<= Tuple(j) .<= (2^l)+1)
    HB.levels[l+1].values[j] = X 
end

function Base.setindex!(NB::NodalBasis{T}, X::T, j::CartesianIndex{D})::T where {T,D}
    @assert all(1 .<= Tuple(j) .<= (2^levels(NB))+1)
    NB.values[j] = X
end

#---------------------------------------------------------------
# Construct basis transformation utilities
#---------------------------------------------------------------

function evaluate(HB::HierarchicalBasis{T,D}, x::CartesianIndex{D})::T where {T,D}
    value = T(0) 
    for l in 0:levels(HB), j in CartesianIndices(ntuple(length->2^l, D))
        value = value + HB[l,j]*basis(l,j,x)
    end
    return value
end

function Base.zeros(NB::Type{NodalBasis{T,D}}, N::Int)::NodalBasis{T,D} where {T,D}
    try 
        Int(log2(N-1)) 
    catch 
        @error "N is not 2^l + 1"
        exit()
    end
    return NodalBasis(zeros(T, ntuple(x->N, D)))
end

function Base.CartesianIndices(NB::NodalBasis{T,D}) where {T,D}
    return CartesianIndices(NB.values)
end

function Base.map(NB::Type{NodalBasis{T,D}}, N::Int, u::Function)::NodalBasis{T,D} where {D, T}
    U = zeros(NodalBasis{T,D}, N)
    for j in CartesianIndices(U)
        U[j] = u(x(Int(log2(N-1)), j)...)
    end
    return U
end
