#-----------------------------------------------------
# SparseGrids | Hierarchical Basis
# Soham 07/2020
#-----------------------------------------------------

using LinearAlgebra
export Level, HierarchicalBasis, NodalBasis
export basis, x, levels, mapme

#---------------------------------------------------------------
# Construct basis function and collocation points
#---------------------------------------------------------------

function basis(l::CartesianIndex{D}, j::CartesianIndex{D}, x::NTuple{D,T})::T where {T,D}
    """
    Evaluate basis functions for level (l) and position (j).
    range(l) := [0, ...]
    range(j) := [0, 2^l]
    """
    l = Tuple(l) .- 1
    j = Tuple(j) .- 1

    @assert all(0 .<= l)
    @assert all(0 .<= j .<= (2 .^l))

    h = 1 ./(2 .^l)

	if all(h.*(j .- 1) .<= x .<= h.*(j .+ 1)) && all(0 .<= x .<= 1)
        value = prod(1 .- abs.((x./h) .- j))
	else
		value = 0
	end
	return value
end

function x(::Type{T}, l::CartesianIndex{D}, j::CartesianIndex{D})::NTuple{D, T} where {D,T}
    l = Tuple(l) .- 1 
    j = Tuple(j) .- 1 

    @assert all(0 .<= l)
    @assert all(0 .<= j .<= 2 .^l)

    h = 1 ./(2 .^l)
    return h.*j
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
    levels::Dict{CartesianIndex{D}, Level{T,D}}
end

struct NodalBasis{T,D}
    values::Array{T,D}
end

#---------------------------------------------------------------
# Construct indexing utilities
#---------------------------------------------------------------

function levels(HB::HierarchicalBasis{T,D})::NTuple{D,Int} where {T,D}
     c = length(HB.levels)
     l = Int((1/2)*(-1 + sqrt(1 + 4c)))
     return (l, l)
end

function levels(NB::NodalBasis{T,D})::NTuple{D,Int} where {T,D}
    return Int.(log2.(size(NB.values) .- 1))
end

function Base.getindex(HB::HierarchicalBasis{T,D}, l::CartesianIndex{D}, j::CartesianIndex{D})::T where {T,D}
    @assert all(1 .<= Tuple(l) .<= levels(HB) .+ 1)
    @assert all(1 .<= Tuple(j) .<= (2 .^Tuple(l)) .+ 1)
    return HB.levels[l].values[j]
end

function Base.setindex!(HB::HierarchicalBasis{T}, X::T, l::CartesianIndex{D}, j::CartesianIndex{D}) where {T,D}
    @assert all(1 .<= Tuple(l) .<= levels(HB) .+ 1)
    @assert all(1 .<= Tuple(j) .<= (2 .^Tuple(l)) .+ 1)
    HB.levels[l].values[j] = X 
end

function Base.getindex(NB::NodalBasis{T}, j::CartesianIndex{D})::T where {T,D}
    @assert all(1 .<= Tuple(j) .<= (2 .^levels(NB)) .+ 1)
    return NB.values[j] 
end

function Base.setindex!(NB::NodalBasis{T}, X::T, j::CartesianIndex{D}) where {T,D}
    @assert all(1 .<= Tuple(j) .<= (2 .^levels(NB)) .+ 1)
    NB.values[j] = X
end

function Base.getindex(NB::NodalBasis{T}, l::CartesianIndex{D}, j::CartesianIndex{D})::T where {T,D}
    @assert all(1 .<= Tuple(l) .<= levels(NB) .+ 1)
    @assert all(1 .<= Tuple(j) .<= (2 .^levels(NB)) .+ 1)
    return NB.values[(2 .^(levels(NB) .- Tuple(l)))*j] 
end

#---------------------------------------------------------------
# Construct basis transformation utilities
# FIXME: Change indexing to remove zeros
# TODO: Construct derivative operators
#---------------------------------------------------------------

function Base.CartesianIndices(NB::NodalBasis{T,D}) where {T,D}
    return CartesianIndices(NB.values)
end

function Base.CartesianIndices(HB::HierarchicalBasis{T,D}) where {T,D}
    return CartesianIndices(levels(HB))
end

function Base.CartesianIndices(HB::HierarchicalBasis{T,D}, l::CartesianIndex{D}) where {T, D}
    return CartesianIndices(2 .^Tuple(l) .+ 1)
end

function Base.zeros(NB::Type{NodalBasis{T,D}}, N::NTuple{D,Int})::NodalBasis{T,D} where {T,D}
    @assert all(mod.(N, 2) .== 1)
    return NodalBasis(zeros(T, N))
end

function evaluate(HB::HierarchicalBasis{T,D}, x::NTuple{D,T})::T where {T,D}
    value = T(0) 
    for l in CartesianIndices(HB), j in CartesianIndices(HB, l)
        value = value + HB[l,j]*basis(l,j,x)
    end
    return value
end

function Base.map(NB::Type{NodalBasis{T,D}}, N::NTuple{D,Int}, u::Function)::NodalBasis{T,D} where {T,D}
    U = zeros(NodalBasis{T,D}, N)
    for j in CartesianIndices(U)
        U[j] = u(x(T, CartesianIndex(levels(U) .+ 1), j)...)
    end
    return U
end

function Base.map(HB::Type{HierarchicalBasis{T,D}}, levels::NTuple{D,Int}, u::Function)::HierarchicalBasis{T,D} where {T,D}
    HB = HierarchicalBasis(Dict{CartesianIndex{D}, Level{T,D}}())
    for l in CartesianIndices(levels), j in CartesianIndices(HB, l)
        HB[l,j] = u(x(T,l,j)...) - evaluate(HB, x(T,l,j))  
    end
    return HB
end

function HierarchicalBasis(NB::NodalBasis{T,D})::HierarchicalBasis{T,D} where {T,D}
    HB = Dict{CartesianIndex{D}, Level{T,D}}()
    for l in CartesianIndices(HB), j in CartesianIndices(HB, l)
        HB[l,j] = NB[l,j] - evaluate(HB, x(T,l,j))  
    end
    return HB
end

function NodalBasis(HB::HierarchicalBasis{T,D})::NodalBasis{T,D} where {T,D}
    NB = zeros(NodalBasis{T,D}, N)
    for j in CartesianIndices(U)
        NB[j] = evaluate(HB, x(T, CartesianIndex(levels(U)), j))
    end
    return NB
end
