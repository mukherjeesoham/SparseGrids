#-----------------------------------------------------
# SparseGrids | LinearInterpolation
# Soham 07/2020
#-----------------------------------------------------

function interpolate(u::NTuple{2, T}, v::NTuple{2, T}, x::T)::T where {T}
    """
    See <https://en.wikipedia.org/wiki/Linear_interpolation> 
    """
    x0, y0 = u
    x1, y1 = v
    @assert  x0 <= x <= x1
    return y0 + (x - x0)*(y - y0)/(x1 - x0)
end
