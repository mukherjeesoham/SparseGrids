#-----------------------------------------------------
# SparseGrids | Hierarchical Basis
# Soham 07/2020
#-----------------------------------------------------

using LinearAlgebra
export Level, HierarchicalBasis, NodalBasis
export NB_2_HB, HB_2_NB, levels

#---------------------------------------------------------------
# Construct basis function and collocation points
# TODO: Generalize to arbitrary-sized domains.
#---------------------------------------------------------------

function basis(l::Int, j::Int, x::T)::T where {T}
    """
    Evaluate basis functions for level (l) and position (j).
    range(l) := [0, ...]
    range(j) := [0, 2^l]
    """
	@assert 0 <= l
	@assert 0 <= j <= (2^l)

    h = 1/(2^l)

	if (h*(j - 1) <= x <= h*(j + 1)) && (0 <= x <= 1)
		value = 1 - abs((x/h) - j)
	else
		value = 0
	end
   
	return value
end

function x(::T, l::Int, j::Int)::T where {T}
	@assert 0 <= l
	@assert 0 <= j <= (2^l)
    h = 1/(2^l)
    return h*j
end

#---------------------------------------------------------------
# Construct types
#---------------------------------------------------------------

struct Level{T}
    values::Array{T,1}
end
    
struct HierarchicalBasis{T}
    levels::Array{Level{T}, 1}
end

struct NodalBasis{T}
    values::Array{T,1}
end

#---------------------------------------------------------------
# Construct indexing utilities
# TODO: Current indexing routine stores a lot of zeros. Fix this. 
#---------------------------------------------------------------

function levels(HB::HierarchicalBasis{T})::Int where {T}
    return length(HB.levels) - 1
end

function levels(NB::NodalBasis{T})::Int where {T}
    return log2(length(NB.values) - 1)
end

function Base.getindex(HB::HierarchicalBasis{T}, l::Int, j::Int)::T where {T}
    @assert 0 <= l <= levels(HB)
	@assert 0 <= j <= (2^l)
    return (HB.levels[l+1]).values[j+1]
end

function Base.getindex(NB::NodalBasis{T}, j::Int)::T where {T}
    @assert 0 <= j <= (2^levels(NB))
    return NB.values[j+1] 
end

function Base.getindex(NB::NodalBasis{T}, l::Int, j::Int)::T where {T}
    @assert 0 <= j <= (2^levels(NB))
    return NB.values[2^(levels(NB)-l)*j + 1] 
end

function Base.setindex!(HB::HierarchicalBasis{T}, X::T, l::Int, j::Int)::T where {T}
    @assert 0 <= l <= levels(HB)
	@assert 0 <= j <= (2^l)
    (HB.levels[l+1]).values[j+1] = X 
end

function Base.setindex!(NB::NodalBasis{T}, X::T, j::Int)::T where {T}
    @assert 0 <= j <= (2^levels(NB))
    NB.values[j+1] = X
end

function Base.lastindex(NB::NodalBasis{T})::Int where {T}
    return length(NB.values) - 1 
end

#---------------------------------------------------------------
# Construct basis transformation utilities
# TODO: To go from a nodal to a modal routine, use a recursive linear
# interpolation routine, to make things more efficient. 
#---------------------------------------------------------------

function evaluate(HB::HierarchicalBasis{T}, x::T)::T where {T}
    value = T(0) 
    for l in 0:levels(HB), j in 0:2^l
        value = value + HB[l,j]*basis(l,j,x)
    end
    return value
end

function Base. zeros(NB::NodalBasis{T})::HierarchicalBasis{T} where {T}
    HB = []
    for l in 0:levels(NB)
        HB = vcat(HB, Level(zeros(T, 2^l + 1)))
    end
    return HierarchicalBasis{T}(HB)
end

function Base. zeros(HB::HierarchicalBasis{T})::NodalBasis{T} where {T}
    return NodalBasis(zeros(T, 2^levels(HB) + 1))
end

function NB_2_HB(NB::NodalBasis{T})::HierarchicalBasis{T} where {T}
    HB = zeros(NB)
    for l in 0:levels(NB), j in 0:2^l
        HB[l,j] = NB[l,j] - evaluate(HB, x(T,l,j))  
    end
    return HB
end

function HB_2_NB(HB::HierarchicalBasis{T})::NodalBasis{T} where {T}
    NB = zeros(HB)
    for j in 0:2^levels(HB)
        NB[j] = evaluate(HB, x(T,levels(HB),j))
    end
    return NB
end
