# ----------------------------------------------------
# SparseGrids.jl
# Soham M 06/2020
# Compute heirarchial basis functions on the 
# interval [0, 1]
# ----------------------------------------------------

export hat

function hat(l::Int, j::Int, x::Float64)::Float64
    @assert !(l < 0)    
    @assert (1 <= j <= 2^l + 1)
    @assert (0 <= x <= 1) 
    hl = 1/(2^l)
    xj = (j-1)*hl 
    if (xj - hl) <= x <= (xj + hl)
        return 1 - abs((x-xj)/hl)
    else
        return 0
    end
end
