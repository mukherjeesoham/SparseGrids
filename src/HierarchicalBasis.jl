function basis(l::Int, j::Int, x::Float64)::Float64
    """
    We set up conventions for hierarchical basis functions here.  Levels go
    from 0 to l, where l = 0 corresponds to the boundary level. j goes from 0
    to 2^l skipping the even points. 
    """
	@assert 0 <= l
	@assert 0 <= j <= (2^l)
    if l > 0
        @assert isodd(j)
    end

	h = 2.0^(-l)
	if (h*(j - 1) <= x <= h*(j + 1)) && (0 <= x <= 1)
		value = 1 - abs((x/h) - j)
	else
		value = 0
	end
   
	return value
end

function x(l::Int, j::Int)::Float64
    h = 2^(-l)
    return h*j
end

struct Level{T}
    values::Array{T,1}
end
    
struct HierarchicalBasis{T}
    levels::Array{Level{T}, 1}
end

struct NodalBasis{T}
    values::Array{T,1}
end

function maxlevel(HB::HierarchicalBasis{T})::Int where {T}
    return length(H.levels) - 1
end

function maxlevel(NB::HierarchicalBasis{T})::Int where {T}
    return log2(length(NB.values) - 1)
end

function Base.getindex(HB::HierarchicalBasis{T}, l::Int, j::Int)::T where {T}
    @assert 0 <= l <= maxlevel(HB)
	@assert 0 <= j <= (2^l)
    return (HB.levels[l+1]).values[div(j+1,2)]
end

function Base.getindex(NB::NodalBasis{T}, j::Int)::T where {T}
    @assert 0 <= j <= (2^maxlevel(NB))
    return NB.values[j+1] 
end

function Base.setindex!(HB::HierarchicalBasis{T}, X::T, l::Int, j::Int)::T where {T}
    @assert 0 <= l <= maxlevel(HB)
	@assert 0 <= j <= (2^l)
    (HB.levels[l+1]).values[div(j+1,2)] = X
end

function Base.setindex!(NB::NodalBasis{T}, X::T, j::Int)::T where {T}
    @assert 0 <= j <= (2^maxlevel(NB))
    NB.values[j+1] = X
end

function evaluate(HB::HierarchicalBasis{T}, x::T)::T where {T}
    value = HB[0, 0]*basis(0, 0, x) + HB[0, 1]*basis(0, 1, x) 
    for l in 1:maxlevel(HB), j in 1:2:2^l
        value = value + HB[l,j]*basis(l,j,x)
    end
    return value
end

function Base. zeros(NB::NodalBasis{T})::HierarchicalBasis{T} where {T}
    HB = [Level([0,0])]
    for l in 1:maxlevel(NB)
        append!(HB, Level(zeros(T, l)))
    end
    return HierarchialBasis(HB)
end

function Base. zeros(HB::HierarchicalBasis{T})::NodalBasis{T} where {T}
    return NodalBasis(zeros(T, maxlevel(HB)^2 + 1))
end


function NB_2_HB(NB::NodalBasis{T})::HierarchicalBasis{T} where {T}
    HB = zeros(NB)
    HB[0,1] = NB[1]
    HB[0,2] = NB[end]
    for l in 1:maxlevel(NB), j in 1:2:2^l
        HB[l,j] = NB[j] - evaluate(HB, x(l,j))  
    end
    return H 
end

function HB_2_NB(NB::HierarchicalBasis{T})::NodalBasis{T} where {T}
    NB = zeros(HB)
    NB[0,1]   = HB[0,1]
    NB[0,end] = HB[0,2]
    for j in 1:2^l - 1
        xj    = (2^(-l))*j
        NB[j] = evaluate(HB, x(l,j))
    end
    return NB
end
